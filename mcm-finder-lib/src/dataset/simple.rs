use std::{
    collections::HashMap,
    fmt::{Display, Formatter},
    fs::File,
    hash::BuildHasher,
    io::{BufReader, Read},
    marker::PhantomData,
    path::Path,
};

use fixedbitset::FixedBitSet;
use miette::NamedSource;
use statrs::function::gamma::ln_gamma;

use super::{Dataset, DefaultState, line_length_tracker, verify_ascii};
use crate::{
    dataset::datacontainer::{DataContainer, DataMap, DataVec},
    mcm_error::MCMError,
};

#[derive(Debug, Clone)]
/// A simple dataset without any caching or optimization.
///
/// # Examples
/// ```
/// use mcm_finder_lib::dataset::{SimpleDataset, datacontainer::DataVec, DefaultState, Dataset};
/// let dataset = SimpleDataset::<DataVec<DefaultState>, DefaultState>::read_from_file(
///     std::path::Path::new("tests/data/SCOTUS_n9_N895_Data.dat"),
/// ).unwrap();
/// assert_eq!(dataset.datapoints(), 895);
/// assert_eq!(dataset.variables(), 9);
/// ```
pub struct SimpleDataset<C: DataContainer<S>, S: BuildHasher + Default> {
    pub(crate) data: C,
    datapoints: usize,
    _build_hasher: PhantomData<S>,
}

pub type VecDataset = SimpleDataset<DataVec<DefaultState>, DefaultState>;
pub type MapDataset = SimpleDataset<DataMap<DefaultState>, DefaultState>;

impl<C: DataContainer<S>, S: BuildHasher + Default> SimpleDataset<C, S> {
    pub fn new(data: C, datapoints: usize) -> SimpleDataset<C, S> {
        SimpleDataset {
            data,
            datapoints,
            _build_hasher: PhantomData,
        }
    }

    /// Returns the count of datapoints matching the given configuration.
    ///
    /// # Examples
    /// ```
    /// use mcm_finder_lib::dataset::{SimpleDataset, datacontainer::DataVec, DefaultState};
    /// use fixedbitset::FixedBitSet;
    /// use std::collections::HashMap;
    ///
    /// let mut bitset1 = FixedBitSet::with_capacity(3);
    /// bitset1.set(0, true);
    ///
    /// let mut bitset2 = FixedBitSet::with_capacity(3);
    /// bitset2.set(0, true);
    /// bitset2.set(1, true);
    ///
    /// let data = DataVec::from(HashMap::from([
    ///     (bitset1.clone(), 10),
    ///     (bitset2.clone(), 20),
    /// ]));
    /// let dataset = SimpleDataset::<DataVec<DefaultState>, DefaultState>::new(data, 2);
    ///
    /// // Matching configuration returns Some(count)
    /// assert_eq!(dataset.get(&bitset1), Some(10));
    ///
    /// // Non-matching configuration returns None
    /// let mut bitset3 = FixedBitSet::with_capacity(3);
    /// bitset3.set(1, true);
    /// assert_eq!(dataset.get(&bitset3), None);
    /// ```
    pub fn get(&self, configuration: &FixedBitSet) -> Option<usize> {
        self.data
            .iter()
            .find_map(|(d, n)| if &d == configuration { Some(n) } else { None })
    }

    /// Returns an iterator over the datapoints in this dataset.
    ///
    /// # Examples
    /// ```
    /// use mcm_finder_lib::dataset::{SimpleDataset, datacontainer::DataVec, DefaultState};
    /// use fixedbitset::FixedBitSet;
    /// use std::collections::HashMap;
    ///
    /// let mut bitset1 = FixedBitSet::with_capacity(3);
    /// bitset1.set(0, true);
    ///
    /// let mut bitset2 = FixedBitSet::with_capacity(3);
    /// bitset2.set(0, true);
    /// bitset2.set(1, true);
    ///
    /// let data = DataVec::from(HashMap::from([
    ///     (bitset1.clone(), 10),
    ///     (bitset2.clone(), 20),
    /// ]));
    /// let dataset = SimpleDataset::<DataVec<DefaultState>, DefaultState>::new(data, 2);
    ///
    /// // Demonstrate iteration over the dataset and collecting results into a vector
    /// let results: Vec<(FixedBitSet, usize)> = dataset.iter().collect();
    ///
    /// // Include assertion checking iteration length equals dataset size
    /// assert_eq!(results.len(), 2);
    /// ```
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (FixedBitSet, usize)> {
        self.data.iter()
    }

    /// Return the histogram of data for the given ICC.
    ///
    /// This method transforms the dataset into a new `SimpleDataset` where each configuration
    /// is projected onto the bitmask provided by `icc`. If multiple configurations from the
    /// original dataset project to the same configuration in the new one, their counts are summed.
    ///
    /// # Examples
    /// ```
    /// use mcm_finder_lib::dataset::{SimpleDataset, datacontainer::DataVec, DefaultState};
    /// use fixedbitset::FixedBitSet;
    /// use std::collections::HashMap;
    ///
    /// let mut bitset1 = FixedBitSet::with_capacity(3);
    /// bitset1.set(0, true); // bits: {0}
    ///
    /// let mut bitset2 = FixedBitSet::with_capacity(3);
    /// bitset2.set(1, true); // bits: {1}
    ///
    /// let mut bitset3 = FixedBitSet::with_capacity(3);
    /// bitset3.set(0, true);
    /// bitset3.set(2, true); // bits: {0, 2}
    ///
    /// let data = DataVec::from(HashMap::from([
    ///     (bitset1.clone(), 10),
    ///     (bitset2.clone(), 20),
    ///     (bitset3.clone(), 10),
    /// ]));
    /// let dataset = SimpleDataset::<DataVec<DefaultState>, DefaultState>::new(data, 2);
    ///
    /// // Transform to ICC with bits {0, 1} set
    /// let mut icc = FixedBitSet::with_capacity(3);
    /// icc.set(0, true);
    /// icc.set(1, true);
    ///
    /// let transformed = dataset.transform_to_icc(&icc);
    ///
    /// // The transformed dataset contains configurations within the ICC.
    /// assert_eq!(transformed.get(&bitset1), Some(20));
    /// assert_eq!(transformed.get(&bitset2), Some(20));
    ///
    /// // Accessing data via iter() and collecting into a vector
    /// let results: Vec<(FixedBitSet, usize)> = transformed.iter().collect();
    /// assert_eq!(results.len(), 2);
    /// ```
    pub fn transform_to_icc(&self, icc: &FixedBitSet) -> SimpleDataset<C, S> {
        SimpleDataset::new(self.data.to_icc(icc), self.datapoints)
    }
}

impl<C: DataContainer<S>, S: BuildHasher + Default> Dataset for SimpleDataset<C, S> {
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
            .iter_counts()
            .map(|k| ln_gamma((k) as f64 + 0.5) - ln_gamma(0.5))
            .sum::<f64>()
    }

    fn read_from_file(path: &Path) -> Result<SimpleDataset<C, S>, MCMError> {
        let mut data = HashMap::with_hasher(S::default());
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

impl<C: DataContainer<S>, S: BuildHasher + Default> Display for SimpleDataset<C, S> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
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
