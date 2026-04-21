use fixedbitset::FixedBitSet;
use miette::NamedSource;
use std::{
    collections::{HashMap, hash_map::IntoIter},
    fs::File,
    io::{BufReader, Read},
    path::Path,
};

use crate::mcm_error::MCMError;

#[derive(Debug)]
pub struct Dataset {
    pub data: HashMap<FixedBitSet, usize>,
    pub datapoints: usize,
}

impl Dataset {
    pub fn new(data: HashMap<FixedBitSet, usize>, datapoints: usize) -> Dataset {
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

        Ok(Dataset::new(data, datapoints))
    }

    pub fn get(&self, configuration: FixedBitSet) -> Option<usize> {
        self.data.get(&configuration).copied()
    }

    pub fn partition(&self, partition: &FixedBitSet) -> Dataset {
        let new_vectors: Vec<(FixedBitSet, usize)> =
            self.iter().map(|(i, &n)| (i & partition, n)).collect();
        // transformed mu vectors are not guaranteed to be unique, so we want to
        // add them together without loosing information
        let mut partitioned_map = HashMap::new();
        for (o, n) in new_vectors {
            partitioned_map
                .entry(o)
                .and_modify(|v| *v += n)
                .or_insert(n);
        }
        Dataset::new(partitioned_map, self.datapoints)
    }

    pub fn iter(&self) -> std::collections::hash_map::Iter<'_, FixedBitSet, usize> {
        self.data.iter()
    }
}

impl IntoIterator for Dataset {
    type Item = (FixedBitSet, usize);
    type IntoIter = IntoIter<FixedBitSet, usize>;

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
