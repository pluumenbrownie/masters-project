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
    path::Path,
    vec::IntoIter,
};

use crate::mcm_error::MCMError;

pub trait Dataset {
    /// Returns the amount of datapoints in this dataset.
    fn datapoints(&self) -> usize;

    /// Returns the amount of variables in each datapoint.
    fn variables(&self) -> usize;

    /// Returns the amount of bins (unique datapoints) in this dataset.
    fn bins(&self) -> usize;

    /// Return the histogram of data for the given ICC.
    fn transform_to_icc(&self, icc: &FixedBitSet) -> VecDataset;

    /// Returns the prevalence of each state of this variable.
    fn state_prevalence(&self, variable: usize) -> Vec<usize>;

    /// Compute the logarithmic evidence for this ICC.
    fn log_e(&self) -> f64;
}

#[derive(Debug, Clone)]
pub struct VecDataset {
    data: Vec<(FixedBitSet, usize)>,
    datapoints: usize,
}

impl VecDataset {
    pub fn new(data: Vec<(FixedBitSet, usize)>, datapoints: usize) -> VecDataset {
        VecDataset { data, datapoints }
    }

    pub fn read_from_file(path: &Path) -> Result<VecDataset, MCMError> {
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

        Ok(VecDataset::new(data.into_iter().collect(), datapoints))
    }

    pub fn get(&self, configuration: &FixedBitSet) -> Option<usize> {
        self.data
            .iter()
            .find_map(|(d, n)| if d == configuration { Some(n) } else { None })
            .copied()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (FixedBitSet, usize)> {
        self.data.iter()
    }
}

impl Dataset for VecDataset {
    fn variables(&self) -> usize {
        self.data.iter().map(|d| d.0.len()).next().unwrap()
    }

    fn datapoints(&self) -> usize {
        self.datapoints
    }

    /// Returns the amount of bins (unique datapoints) in this dataset.
    fn bins(&self) -> usize {
        self.data.len()
    }

    /// Return the histogram of data for the given ICC.
    fn transform_to_icc(&self, partition: &FixedBitSet) -> VecDataset {
        let new_vectors = self.iter().map(|(i, n)| (i & partition, *n));
        // transformed vectors are not guaranteed to be unique, so we want to
        // add them together without loosing information
        let mut partitioned_map: HashMap<FixedBitSet, usize> = HashMap::new();
        for (o, n) in new_vectors {
            partitioned_map
                .entry(o)
                .and_modify(|v| *v += n)
                .or_insert(n);
        }
        VecDataset::new(partitioned_map.into_iter().collect(), self.datapoints)
    }

    fn state_prevalence(&self, variable: usize) -> Vec<usize> {
        let mut mask = FixedBitSet::with_capacity_and_blocks(self.variables(), [0]);
        mask.set(variable, true);
        let partition: VecDataset = self.transform_to_icc(&mask);
        partition.iter().map(|(_, c)| *c).collect()
    }

    fn log_e(&self) -> f64 {
        self.iter()
            .map(|(_, k)| ln_gamma((*k) as f64 + 0.5) - ln_gamma(0.5))
            .sum::<f64>()
    }
}

impl IntoIterator for VecDataset {
    type Item = (FixedBitSet, usize);
    type IntoIter = IntoIter<(FixedBitSet, usize)>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl Display for VecDataset {
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
pub struct EndsCachedDataset {
    pub data: Vec<VecDataset>,
    variables: usize,
    datapoints: usize,
    bins: usize,
}

impl Dataset for EndsCachedDataset {
    fn variables(&self) -> usize {
        self.variables
    }

    fn datapoints(&self) -> usize {
        self.datapoints
    }

    /// Returns the amount of bins (unique datapoints) in this dataset.
    fn bins(&self) -> usize {
        self.bins
    }

    /// Return the histogram of data for the given ICC.
    fn transform_to_icc(&self, partition: &FixedBitSet) -> VecDataset {
        let ends = Ends::from_icc(partition);
        let location = get_ends_cache_location(ends.as_ref().cloned(), self.variables);
        if location > self.data.len() {
            dbg!(ends);
            println!("{partition}");
            dbg!(self.variables);
        }
        self.data[location].transform_to_icc(partition)
    }

    fn state_prevalence(&self, variable: usize) -> Vec<usize> {
        let mut mask = FixedBitSet::with_capacity_and_blocks(self.variables(), [0]);
        mask.set(variable, true);
        let partition: VecDataset = self.transform_to_icc(&mask);
        partition.iter().map(|(_, c)| *c).collect()
    }

    fn log_e(&self) -> f64 {
        todo!()
    }
}

impl EndsCachedDataset {
    pub fn read_from_file(path: &Path) -> Result<EndsCachedDataset, MCMError> {
        let base = VecDataset::read_from_file(path)?;
        // println!("Base length: {}", base.data.len());
        let mut base_ref_index = 0usize;
        let variables = base.variables();
        let mut data = Vec::with_capacity((variables * (variables + 1)) / 2);

        let mut icc = FixedBitSet::with_capacity_and_blocks(variables, [0]);
        for start in 0..variables {
            icc.set_range(0..start, false);
            icc.set_range(start..variables, true);

            data.push(
                data.get(base_ref_index)
                    .unwrap_or(&base)
                    .transform_to_icc(&icc),
            );
            base_ref_index = data.len() - 1;
            println!(
                "{}: {icc} - {}",
                get_ends_cache_location(Some(Ends::new(start, variables)), variables),
                data[base_ref_index].data.len()
            );
            let mut sub_base_index = base_ref_index;

            for end in (start + 1..variables).rev() {
                icc.set(end, false);
                data.push(data[sub_base_index].transform_to_icc(&icc));

                sub_base_index = data.len() - 1;
                println!(
                    "{}: {icc} - {}",
                    get_ends_cache_location(Some(Ends::new(start, end)), variables),
                    data[sub_base_index].data.len()
                );
            }
        }
        icc.set_range(.., false);
        data.push(data.last().unwrap().transform_to_icc(&icc));
        println!("{}: {icc}", get_ends_cache_location(None, variables));

        Ok(EndsCachedDataset {
            datapoints: data[0].datapoints(),
            bins: data[0].bins(),
            variables: data[0].variables(),
            data,
        })
    }

    pub fn get(&self, configuration: &FixedBitSet) -> Option<usize> {
        self.data[0]
            .iter()
            .find_map(|(d, n)| if d == configuration { Some(n) } else { None })
            .copied()
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

// impl Dataset for EndsCachedDataset {

// }
