use std::{
    fmt::{Display, Formatter},
    hash::BuildHasher,
    marker::PhantomData,
    path::Path,
};

use fixedbitset::FixedBitSet;
use statrs::function::gamma::ln_gamma;

use super::{Dataset, DefaultState, simple::SimpleDataset};
use crate::{
    dataset::{
        datacontainer::{DataContainer, DataVec},
        ends::{Ends, get_ends_cache_location},
    },
    mcm_error::MCMError,
};

#[derive(Debug, Clone)]
/// A simple dataset without any caching or optimization.
///
/// # Examples
/// ```
/// use mcm_finder_lib::dataset::{SimpleDataset, DataVec, DefaultState, Dataset};
/// let dataset = SimpleDataset::<DataVec<DefaultState>, DefaultState>::read_from_file(
///     std::path::Path::new("tests/data/SCOTUS_n9_N895_Data.dat"),
/// ).unwrap();
/// assert_eq!(dataset.datapoints(), 895);
/// assert_eq!(dataset.variables(), 9);
/// ```
pub struct EndsCachedDataset<C: DataContainer<S>, S: BuildHasher + Default> {
    pub data: Vec<Option<C>>,
    pub(crate) variables: usize,
    pub(crate) datapoints: usize,
    pub(crate) _build_hasher: PhantomData<S>,
}

pub type EndsCachedVecDataset = EndsCachedDataset<DataVec<DefaultState>, DefaultState>;

impl<C: DataContainer<S>, S: BuildHasher + Default> Dataset for EndsCachedDataset<C, S> {
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

    fn read_from_file(path: &Path) -> Result<EndsCachedDataset<C, S>, MCMError> {
        let base = SimpleDataset::read_from_file(path)?;
        // println!("Base length: {}", base.data.len());
        let variables = base.variables();
        let data: Vec<Option<C>> = Vec::with_capacity((variables * (variables + 1)) / 2);

        let mut output = EndsCachedDataset {
            datapoints: base.datapoints(),
            variables: base.variables(),
            data,
            _build_hasher: PhantomData,
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

impl<C: DataContainer<S>, S: BuildHasher + Default> Display for EndsCachedDataset<C, S> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
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

impl<C: DataContainer<S>, S: BuildHasher + Default> EndsCachedDataset<C, S> {
    /// Get the amount of counts for this ICC.
    ///
    /// Returns `None` when the dataset is empty.
    ///
    /// # Examples
    /// ```
    /// use mcm_finder_lib::dataset::{EndsCachedVecDataset, DefaultState, Dataset, DataContainer, DataVec};
    /// use fixedbitset::FixedBitSet;
    ///
    /// let dataset = EndsCachedVecDataset::read_from_file(
    ///     std::path::Path::new("tests/data/SCOTUS_n9_N895_Data.dat")
    /// ).unwrap();
    ///
    /// let mut icc = FixedBitSet::with_capacity(9);
    /// icc.set(2, true);
    ///
    /// let selection = dataset.get_icc(&icc).unwrap();
    /// assert_eq!(selection, 1);
    /// ```
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
    ///
    /// # Examples
    /// ```
    /// use mcm_finder_lib::dataset::{EndsCachedVecDataset, DefaultState, Dataset, DataContainer, DataVec};
    /// use fixedbitset::FixedBitSet;
    ///
    /// let dataset = EndsCachedVecDataset::read_from_file(
    ///     std::path::Path::new("tests/data/SCOTUS_n9_N895_Data.dat")
    /// ).unwrap();
    ///
    /// let selection = dataset.get(27).unwrap();
    /// assert_eq!(selection.bins(), 8);
    ///
    /// let none_selection = dataset.get(0).unwrap();
    /// assert_eq!(none_selection.bins(), 128);
    /// ```
    pub fn get(&self, index: usize) -> Option<&C> {
        if self.data.is_empty() {
            return None;
        }
        self.data[..=index].iter().rev().find_map(|d| d.as_ref())
    }

    /// Turns this EndsCachedDataset into a SimpleDataset, transformed for the given ICC.
    ///
    /// # Examples
    /// ```
    /// use mcm_finder_lib::dataset::{EndsCachedVecDataset, DefaultState, Dataset};
    /// use fixedbitset::FixedBitSet;
    ///
    /// let dataset = EndsCachedVecDataset::read_from_file(
    ///     std::path::Path::new("tests/data/SCOTUS_n9_N895_Data.dat")
    /// ).unwrap();
    ///
    /// let mut icc = FixedBitSet::with_capacity(10);
    /// icc.set(0, true);
    /// icc.set(1, true);
    ///
    /// let transformed = dataset.transform_to_icc(&icc);
    ///
    /// for (data, count) in transformed.iter() {
    ///     assert!(count > 0);
    /// }
    /// ```
    pub fn transform_to_icc(&self, partition: &FixedBitSet) -> SimpleDataset<C, S> {
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
