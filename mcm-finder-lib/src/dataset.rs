use fixedbitset::FixedBitSet;
use miette::NamedSource;
use statrs::function::gamma::ln_gamma;
use std::{
    collections::HashMap,
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
    fn transform_to_icc(&self, icc: &FixedBitSet) -> impl LogE;
}

pub trait LogE {
    /// Compute the logarithmic evidence for this dataset.
    fn log_e(&self) -> f64;
}

#[derive(Debug)]
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

    // /// Return the histogram of data for the given ICC.
    // fn transform_to_icc(&self, partition: &FixedBitSet) -> VecDataset {
    //     let new_vectors = self.iter().map(|(i, n)| (i & partition, *n));
    //     // transformed vectors are not guaranteed to be unique, so we want to
    //     // add them together without loosing information
    //     let mut partitioned_map: Vec<(FixedBitSet, usize)> = Vec::new();
    //     for (o, n) in new_vectors {
    //         let exists = partitioned_map.binary_search_by(|(b, _)| b.cmp(&o));
    //         match exists {
    //             Ok(i) => partitioned_map[i].1 += n,
    //             Err(i) => {
    //                 if i == partitioned_map.len() {
    //                     partitioned_map.push((o, n));
    //                 } else {
    //                     partitioned_map.insert(i, (o, n));
    //                 }
    //             }
    //         }
    //     }
    //     VecDataset::new(partitioned_map.into_iter().collect(), self.datapoints)
    // }

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
}

impl LogE for VecDataset {
    fn log_e(&self) -> f64 {
        self.iter()
            .map(|(_, k)| ln_gamma(*k as f64 + 0.5) - ln_gamma(0.5))
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

// #[derive(Debug, Clone)]
// pub struct TreeDataset {
//     datapoints: usize,
//     branches: Vec<Option<TreeDataset>>,
// }

// impl TreeDataset {
//     pub fn read_from_file(path: &Path) -> Result<TreeDataset, MCMError> {
//         let mut tree: Option<TreeDataset> = None;
//         let mut datapoints = 0usize;
//         let filename = path.file_name().unwrap().to_str().unwrap().to_owned();
//         let mut buf_reader = BufReader::new(File::open(path)?);

//         // reading the entire file to a string
//         let file = {
//             let mut file = String::new();
//             buf_reader.read_to_string(&mut file)?;
//             file
//         };

//         let mut line_length = 0usize;
//         let mut bool_array: Vec<bool> = vec![];

//         for (nr, byte) in file.bytes().enumerate() {
//             // validate character is valid ascii
//             verify_ascii(&filename, &file, nr, byte)?;
//             // count ones and zeroes
//             match byte {
//                 b'0' => bool_array.push(false),
//                 b'1' => bool_array.push(true),
//                 b'\r' | b'\n' if !bool_array.is_empty() => {
//                     // check the line length
//                     line_length_tracker(&filename, &file, &mut line_length, &bool_array, nr)?;

//                     // add the datapoints to the hashmap
//                     let mut bitvec = FixedBitSet::with_capacity(bool_array.len());
//                     for (nr, bit) in bool_array.drain(..).enumerate() {
//                         bitvec.set(nr, bit);
//                     }
//                     tree.get_or_insert(TreeDataset::new(line_length))
//                         .add_bitvector(bitvec);
//                     // data.entry(bitvec).and_modify(|i| *i += 1).or_insert(1usize);
//                     debug_assert!(bool_array.is_empty());

//                     datapoints += 1;
//                 }
//                 b'\r' | b'\n' => {}
//                 // wrong character case
//                 _ => Err(MCMError::BadCharacter {
//                     src: NamedSource::new(&filename, file.clone()),
//                     bad_line: nr.into(),
//                 })?,
//             }
//         }

//         Ok(tree.unwrap())
//     }

//     /// Add a datapoint to this dataset.
//     pub fn add_bitvector(&mut self, datapoint: FixedBitSet) {
//         self.add(datapoint.into_ones().rev().collect());
//     }

//     /// Add a datapoint to this dataset in vector form. Vector must be sorted in
//     /// decrementing order.
//     pub fn add(&mut self, datapoint: Vec<usize>) {
//         let mut datapoint = datapoint;
//         self.datapoints += 1;
//         if let Some(nr) = datapoint.pop() {
//             if self.branches[nr].is_none() {
//                 self.branches[nr] = Some(TreeDataset::new(self.variables() - nr - 1));
//             }
//             self.branches[nr]
//                 .as_mut()
//                 .unwrap()
//                 .add(datapoint.iter_mut().map(|p| *p - nr - 1).collect());
//             // match self.branches[nr] {
//             //     None => self.branches[nr] = Some(TreeDataset::one(self.variables() - nr)),
//             //     Some(branch) => branch.add(datapoint),
//             // }
//         }
//     }

//     pub fn new(variables: usize) -> TreeDataset {
//         TreeDataset {
//             datapoints: 0,
//             branches: vec![None; variables],
//         }
//     }
// }

// impl Dataset for TreeDataset {
//     fn datapoints(&self) -> usize {
//         self.datapoints
//     }

//     fn bins(&self) -> usize {
//         if self.branches.is_empty() {
//             1
//         } else {
//             1 + self
//                 .branches
//                 .iter()
//                 .map(|b| match b {
//                     None => 0,
//                     Some(t) => t.bins(),
//                 })
//                 .sum::<usize>()
//         }
//     }

//     fn variables(&self) -> usize {
//         self.branches.len()
//     }

//     fn transform_to_icc(&self, icc: &FixedBitSet) -> impl LogE {
//         todo!();
//         VecDataset::new(vec![], 0)
//     }
// }

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
