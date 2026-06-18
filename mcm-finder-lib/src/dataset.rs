//! This module contains the different kinds of datasets and relevant traits in
//! this library.

use fixedbitset::FixedBitSet;
use miette::NamedSource;
use statrs::function::gamma::ln_gamma;
use std::{
    collections::HashMap,
    fmt::{self, Display},
    fs::File,
    io::{BufReader, Read},
    marker::PhantomData,
    path::Path,
    vec::IntoIter,
};

use crate::mcm_error::MCMError;

// trait BitVector {}

// impl BitVector for FixedBitSet {}

pub trait DataContainter: From<HashMap<FixedBitSet, usize>> {
    fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)>;

    fn to_icc(&self, icc: &FixedBitSet) -> Self {
        let mut new_icc_map: HashMap<FixedBitSet, usize> = HashMap::new();
        for (o, n) in self.iter().map(|(i, n)| (i & icc, n)) {
            new_icc_map.entry(o).and_modify(|v| *v += n).or_insert(n);
        }
        Self::from(new_icc_map)
    }

    fn bins(&self) -> usize {
        self.iter().len()
    }
}

#[derive(Debug, Clone)]
pub struct DataVec {
    data: Vec<(FixedBitSet, usize)>,
}

impl From<HashMap<FixedBitSet, usize>> for DataVec {
    fn from(value: HashMap<FixedBitSet, usize>) -> Self {
        Self {
            data: value.into_iter().collect(),
        }
    }
}

impl DataContainter for DataVec {
    fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)> {
        self.data.iter().map(|(v, n)| (v, *n))
    }
}

#[derive(Debug, Clone)]
pub struct DataMap {
    data: HashMap<FixedBitSet, usize>,
}

impl From<HashMap<FixedBitSet, usize>> for DataMap {
    fn from(value: HashMap<FixedBitSet, usize>) -> Self {
        Self { data: value.into() }
    }
}

impl DataContainter for DataMap {
    fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)> {
        self.data.iter().map(|(v, n)| (v, *n))
    }
}

pub trait Dataset {
    /// Returns the amount of datapoints in this dataset.
    fn datapoints(&self) -> usize;

    /// Returns the amount of variables in each datapoint.
    fn variables(&self) -> usize;

    /// Returns the prevalence of each state of this variable.
    fn state_prevalence(&self, variable: usize) -> Vec<usize>;

    /// Compute the logarithmic evidence for the given ICC.
    fn log_e(&self, icc: &FixedBitSet) -> f64;

    fn read_from_file(path: &Path) -> Result<Self, MCMError>
    where
        Self: Sized;
}

#[derive(Debug, Clone)]
pub struct SimpleDataset<C: DataContainter> {
    data: C,
    datapoints: usize,
}

pub type VecDataset = SimpleDataset<DataVec>;
pub type MapDataset = SimpleDataset<DataMap>;

impl<C: DataContainter> SimpleDataset<C> {
    fn new(data: C, datapoints: usize) -> SimpleDataset<C> {
        SimpleDataset { data, datapoints }
    }

    pub fn get(&self, configuration: &FixedBitSet) -> Option<usize> {
        self.data
            .iter()
            .find_map(|(d, n)| if d == configuration { Some(n) } else { None })
        // .copied()
    }

    pub fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)> {
        self.data.iter()
    }

    /// Return the histogram of data for the given ICC.
    pub fn transform_to_icc(&self, icc: &FixedBitSet) -> SimpleDataset<C> {
        SimpleDataset::new(self.data.to_icc(icc), self.datapoints)
    }
}

impl<C: DataContainter> Dataset for SimpleDataset<C> {
    fn variables(&self) -> usize {
        self.data.iter().map(|d| d.0.len()).next().unwrap()
    }

    fn datapoints(&self) -> usize {
        self.datapoints
    }

    fn state_prevalence(&self, variable: usize) -> Vec<usize> {
        let mut mask = FixedBitSet::with_capacity_and_blocks(self.variables(), [0]);
        mask.set(variable, true);
        let partition = self.data.to_icc(&mask);
        partition.iter().map(|(_, c)| c).collect()
    }

    fn log_e(&self, icc: &FixedBitSet) -> f64 {
        self.data
            .to_icc(icc)
            .iter()
            .map(|(_, k)| ln_gamma((k) as f64 + 0.5) - ln_gamma(0.5))
            .sum::<f64>()
    }

    fn read_from_file(path: &Path) -> Result<SimpleDataset<C>, MCMError> {
        let mut data = HashMap::new();
        let mut datapoints = 0usize;
        let filename = path.file_name().unwrap().to_str().unwrap().to_owned();
        let mut buf_reader = BufReader::new(File::open(path)?);

        // reading the entire file to a string
        let file = {
            let mut file = String::new();
            buf_reader.read_to_string(&mut file)?;
            file
        };

        let mut line_length = 0usize;
        let mut bool_array: Vec<bool> = vec![];

        for (nr, byte) in file.bytes().enumerate() {
            // validate character is valid ascii
            verify_ascii(&filename, &file, nr, byte)?;
            // count ones and zeroes
            match byte {
                b'0' => bool_array.push(false),
                b'1' => bool_array.push(true),
                b'\r' | b'\n' if !bool_array.is_empty() => {
                    // check the line length
                    line_length_tracker(&filename, &file, &mut line_length, &bool_array, nr)?;

                    // add the datapoints to the hashmap
                    let mut bitvec = FixedBitSet::with_capacity(bool_array.len());
                    for (nr, bit) in bool_array.drain(..).enumerate() {
                        bitvec.set(nr, bit);
                    }
                    data.entry(bitvec).and_modify(|i| *i += 1).or_insert(1usize);
                    debug_assert!(bool_array.is_empty());

                    datapoints += 1;
                }
                b'\r' | b'\n' => {}
                // wrong character case
                _ => Err(MCMError::BadCharacter {
                    src: NamedSource::new(&filename, file.clone()),
                    bad_line: nr.into(),
                })?,
            }
        }

        Ok(SimpleDataset::new(C::from(data), datapoints))
    }
}

// impl<C: DataContainter> IntoIterator for VecDataset<C> {
//     type Item = (FixedBitSet, usize);
//     type IntoIter = IntoIter<(FixedBitSet, usize)>;

//     fn into_iter(self) -> Self::IntoIter {
//         self.data.into_iter()
//     }
// }

impl<C: DataContainter> Display for SimpleDataset<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "VecDataset: {{")?;
        writeln!(f, "   data: [")?;
        for (v, n) in self.data.iter() {
            writeln!(f, "       ({}, {}),", v, n)?;
        }
        writeln!(f, "   ],")?;
        writeln!(f, "   datapoints: {}", self.datapoints)?;
        writeln!(f, "}}")
    }
}

pub(crate) fn line_length_tracker(
    filename: &str,
    file: &str,
    line_length: &mut usize,
    bool_array: &[bool],
    nr: usize,
) -> Result<(), MCMError> {
    // set the line length if it hasn't been set yet
    if *line_length == 0 {
        *line_length = bool_array.len();
    } else if bool_array.len() != *line_length {
        Err(MCMError::BadLength {
            src: NamedSource::new(filename, file.to_owned()),
            bad_line: (nr - bool_array.len(), bool_array.len()).into(),
        })?
    };
    Ok(())
}

pub(crate) fn verify_ascii(
    filename: &str,
    file: &str,
    char_nr: usize,
    byte: u8,
) -> Result<(), MCMError> {
    if !byte.is_ascii() {
        Err(MCMError::BadCharacter {
            src: NamedSource::new(filename, file.to_owned()),
            bad_line: char_nr.into(),
        })?
    };
    Ok(())
}

#[derive(Debug, Clone)]
pub struct EndsCachedDataset<C: DataContainter> {
    pub data: Vec<Option<C>>,
    variables: usize,
    datapoints: usize,
}

pub type EndsCachedVecDataset = EndsCachedDataset<DataVec>;

impl<C: DataContainter> Dataset for EndsCachedDataset<C> {
    fn variables(&self) -> usize {
        self.variables
    }

    fn datapoints(&self) -> usize {
        self.datapoints
    }

    fn state_prevalence(&self, variable: usize) -> Vec<usize> {
        let mut mask = FixedBitSet::with_capacity_and_blocks(self.variables(), [0]);
        mask.set(variable, true);
        let partition: C = self
            .get(get_ends_cache_location(
                Ends::from_icc(&mask),
                self.variables,
            ))
            .unwrap()
            .to_icc(&mask);
        partition.iter().map(|(_, c)| c).collect()
    }

    fn log_e(&self, icc: &FixedBitSet) -> f64 {
        self.get(get_ends_cache_location(Ends::from_icc(icc), self.variables))
            .unwrap()
            .to_icc(icc)
            .iter()
            .map(|(_, k)| ln_gamma((k) as f64 + 0.5) - ln_gamma(0.5))
            .sum::<f64>()
    }

    fn read_from_file(path: &Path) -> Result<EndsCachedDataset<C>, MCMError> {
        let base = SimpleDataset::read_from_file(path)?;
        // println!("Base length: {}", base.data.len());
        let variables = base.variables();
        let data: Vec<Option<C>> = Vec::with_capacity((variables * (variables + 1)) / 2);

        let mut output = EndsCachedDataset {
            datapoints: base.datapoints(),
            variables: base.variables(),
            data,
        };

        let mut base_ref_index = 0usize;

        let mut icc = FixedBitSet::with_capacity_and_blocks(variables, [0]);
        for start in 0..variables {
            icc.set_range(0..start, false);
            icc.set_range(start..variables, true);

            output.data.push(Some(
                output
                    .get(base_ref_index)
                    .unwrap_or(&base.data)
                    .to_icc(&icc),
            ));
            base_ref_index = output.data.len() - 1;
            // println!(
            //     "{}: {icc} - {}",
            //     get_ends_cache_location(Some(Ends::new(start, variables)), variables),
            //     data[base_ref_index].data.len()
            // );
            let mut sub_base_index = base_ref_index;

            for end in (start + 1..variables).rev() {
                icc.set(end, false);
                if icc.count_ones(..) < variables / 2 {
                    output
                        .data
                        .push(Some(output.get(sub_base_index).unwrap().to_icc(&icc)));
                } else {
                    output.data.push(None);
                }

                sub_base_index = output.data.len() - 1;
                // println!(
                //     "{}: {icc} - {}",
                //     get_ends_cache_location(Some(Ends::new(start, end)), variables),
                //     data[sub_base_index].data.len()
                // );
            }
        }
        icc.set_range(.., false);
        output.data.push(Some(
            output.data.last().unwrap().as_ref().unwrap().to_icc(&icc),
        ));
        // println!("{}: {icc}", get_ends_cache_location(None, variables));

        Ok(output)
    }
}

impl<C: DataContainter> Display for EndsCachedDataset<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut icc = FixedBitSet::with_capacity_and_blocks(self.variables(), [0]);
        let mut data_iter = self.data.iter();
        let variables = self.variables();
        for start in 0..variables {
            icc.set_range(0..start, false);
            icc.set_range(start..variables, true);

            match data_iter.next().unwrap() {
                Some(vec_dataset) => writeln!(
                    f,
                    "{}: {icc} - {}",
                    get_ends_cache_location(Some(Ends::new(start, variables)), variables),
                    vec_dataset.iter().len()
                )?,
                None => writeln!(f, "None")?,
            };

            for end in (start + 1..variables).rev() {
                icc.set(end, false);
                match data_iter.next().unwrap() {
                    Some(vec_dataset) => writeln!(
                        f,
                        "{}: {icc} - {}",
                        get_ends_cache_location(Some(Ends::new(start, end)), variables),
                        vec_dataset.iter().len()
                    )?,
                    None => writeln!(f, "None")?,
                };
            }
        }
        icc.set_range(.., false);
        writeln!(
            f,
            "{}: {icc} - {}",
            get_ends_cache_location(None, variables),
            data_iter.next().unwrap().as_ref().unwrap().iter().len()
        )?;
        Ok(())
    }
}

impl<C: DataContainter> EndsCachedDataset<C> {
    pub fn get_icc(&self, configuration: &FixedBitSet) -> Option<usize> {
        self.data[0]
            .as_ref()
            .unwrap()
            .iter()
            .find_map(|(d, n)| if d == configuration { Some(n) } else { None })
    }

    /// Get the first best dataset for this `get_ends_cache_location` index.
    ///
    /// Returns `None` when the dataset is empty.
    pub fn get(&self, index: usize) -> Option<&C> {
        if self.data.is_empty() {
            return None;
        }
        self.data[..=index].iter().rev().find_map(|d| d.as_ref())
    }

    /// Turns this EndsCachedDataset into a SimpleDataset, transformed for the given ICC.
    pub fn transform_to_icc(&self, partition: &FixedBitSet) -> SimpleDataset<C> {
        let ends = Ends::from_icc(partition);
        let location = get_ends_cache_location(ends.as_ref().cloned(), self.variables);
        if location > self.data.len() {
            dbg!(ends);
            println!("{partition}");
            dbg!(self.variables);
        }
        SimpleDataset::new(
            self.get(location).unwrap().to_icc(partition),
            self.datapoints,
        )
    }
}

/// The range of an ICC which contains included variables.
///
/// - `start` is the index of the first included variable.
///
/// - `end` is the index of the first excluded variable after the last included variable.
#[derive(Debug, Clone)]
struct Ends {
    start: usize,
    end: usize,
}

impl Ends {
    fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    fn from_icc(partition: &FixedBitSet) -> Option<Self> {
        partition.minimum().map(|start| Ends {
            start,
            end: partition.maximum().unwrap() + 1,
        })
    }
}

const fn get_ends_cache_location(ends: Option<Ends>, variables: usize) -> usize {
    match ends {
        None => (variables * (variables + 1)) / 2,
        Some(ends) => {
            if ends.start == 0 && ends.end == variables + 1 {
                0
            } else {
                let start_offset = variables - ends.start;
                (variables.pow(2) + variables - (start_offset.pow(2) + start_offset)) / 2
                    + variables
                    - ends.end
            }
        }
    }
}
