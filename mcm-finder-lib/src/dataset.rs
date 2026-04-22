use fixedbitset::FixedBitSet;
use miette::NamedSource;
use statrs::function::gamma::ln_gamma;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Read},
    num::{NonZeroU8, NonZeroUsize},
    path::Path,
    vec::IntoIter,
};

use crate::mcm_error::MCMError;

struct Branch {
    value: NonZeroUsize,
    next: Vec<Amalia>,
}

enum Amalia {
    None,
    Value(NonZeroU8),
    Branch(Branch),
}

fn amalia() -> Amalia {
    Amalia::None
}

#[derive(Debug)]
pub struct Dataset {
    data: Vec<(FixedBitSet, usize)>,
    datapoints: usize,
}

impl Dataset {
    pub fn new(data: Vec<(FixedBitSet, usize)>, datapoints: usize) -> Dataset {
        Dataset { data, datapoints }
    }

    pub fn read_from_file(path: &Path) -> Result<Dataset, MCMError> {
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

        Ok(Dataset::new(data.into_iter().collect(), datapoints))
    }

    pub fn variables(&self) -> usize {
        self.data.iter().map(|d| d.0.len()).next().unwrap()
    }

    pub fn datapoints(&self) -> usize {
        self.datapoints
    }

    /// Returns the amount of bins (unique datapoints) in this dataset.
    pub fn bins(&self) -> usize {
        self.data.len()
    }

    pub fn get(&self, configuration: &FixedBitSet) -> Option<usize> {
        self.data
            .iter()
            .find_map(|(d, n)| if d == configuration { Some(n) } else { None })
            .copied()
    }

    /// Return the histogram of data for the given ICC.
    pub fn partition(&self, partition: &FixedBitSet) -> Dataset {
        let new_vectors = self.iter().map(|(i, n)| (i & partition, *n));
        // transformed vectors are not guaranteed to be unique, so we want to
        // add them together without loosing information
        let mut partitioned_map: Vec<(FixedBitSet, usize)> = Vec::new();
        for (o, n) in new_vectors {
            let exists = partitioned_map.binary_search_by(|(b, _)| b.cmp(&o));
            match exists {
                Ok(i) => partitioned_map[i].1 += n,
                Err(i) => {
                    if i == partitioned_map.len() {
                        partitioned_map.push((o, n));
                    } else {
                        partitioned_map.insert(i, (o, n));
                    }
                }
            }
        }
        Dataset::new(partitioned_map.into_iter().collect(), self.datapoints)
    }

    pub fn log_e(&self) -> f64 {
        self.iter()
            .map(|(_, k)| ln_gamma(*k as f64 + 0.5) - ln_gamma(0.5))
            .sum::<f64>()
    }

    // pub fn partition_log_e(&self, icc: &FixedBitSet) -> f64 {

    // }

    // pub fn partition(&self, partition: &FixedBitSet) -> Dataset {
    //     let new_vectors = self.iter().map(|(i, n)| (i & partition, *n));
    //     // transformed vectors are not guaranteed to be unique, so we want to
    //     // add them together without loosing information
    //     let mut partitioned_map = Vec::new();
    //     for (o, n) in new_vectors {
    //         let point = partitioned_map.partition_point(|(x, _)| *x < o);
    //         if point == partitioned_map.len() {
    //             partitioned_map.push((o, n));
    //         } else if partitioned_map[point].0 == o {
    //             partitioned_map[point].1 += n;
    //         } else {
    //             partitioned_map.insert(point, (o, n));
    //         }
    //     }
    //     Dataset::new(partitioned_map.into_iter().collect(), self.datapoints)
    // }

    pub fn iter(&self) -> std::slice::Iter<'_, (FixedBitSet, usize)> {
        self.data.iter()
    }
}

impl IntoIterator for Dataset {
    type Item = (FixedBitSet, usize);
    type IntoIter = IntoIter<(FixedBitSet, usize)>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
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
