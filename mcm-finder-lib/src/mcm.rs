//! This module contains the main MinimallyComplexModel type.

use annolog::CollectorEvent;
use dashmap::DashMap;
use fixedbitset::FixedBitSet;
use miette::NamedSource;
use rand::{
    RngExt, SeedableRng,
    distr::Distribution,
    seq::{IndexedRandom, IteratorRandom},
};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use serde::Serialize;
use statrs::{distribution::Binomial, function::gamma::ln_gamma};
use std::{
    collections::{HashMap, VecDeque},
    f64::consts::{LN_2, PI},
    fmt::Display,
    fs::File,
    io::{BufReader, Read},
    mem,
    num::{NonZeroU32, NonZeroUsize},
    ops::Deref,
    path::Path,
    sync::Arc,
};

use crate::{
    dataset::{Dataset, line_length_tracker, verify_ascii},
    logger::SolverEvent,
    mcm::icc::IndependentCompleteComponent,
    mcm_error::MCMError,
};

pub mod icc;

#[derive(Debug, Clone, Copy, Serialize)]
pub enum MutationType {
    Split,
    Merge,
    Swap,
}

impl MutationType {
    /// Returns a random arm of this enum. `merge_prob` gives the probability of
    /// returning `MutationType::Merge`.
    fn rand(rng: &mut rand::rngs::ThreadRng, merge_prob: [f64; 3]) -> MutationType {
        [
            (MutationType::Merge, merge_prob[0]),
            (MutationType::Split, merge_prob[1]),
            (MutationType::Swap, merge_prob[2]),
        ]
        .choose_weighted(rng, |item| item.1)
        .map(|t| t.0)
        .unwrap()
    }
}
#[derive(Debug, Clone, Copy)]
pub enum MCMData {
    Mutation(MutationEvent),
    LogE(f64),
}

impl From<MCMData> for CollectorEvent<SolverEvent> {
    fn from(value: MCMData) -> Self {
        CollectorEvent::Data(SolverEvent::Mcm(value))
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct MutationEvent {
    pub mut_type: MutationType,
    pub accepted: bool,
    pub temperature: f64,
}

/// Calculates the geometric complexity of an Independent Complete Component.
///
/// The geometric complexity of a Minimally Complex Model is equal to the
/// sum of the geometric complexity of its ICCs, with the complexity of an ICC being:
///
/// > 2ʳᵃ⁻¹ log π - log Γ(2ʳᵃ⁻¹)
///
/// with rₐ the modeled spin variables.
///
/// # Example
/// ```
/// # use mcm_finder_lib::mcm::geometric_complexity_icc;
/// # use std::num::NonZeroU32;
/// let spin_variables = NonZeroU32::new(9).unwrap();
/// assert_eq!(geometric_complexity_icc(spin_variables), -868.6612503409542);
/// ```
pub fn geometric_complexity_icc(spin_variables: NonZeroU32) -> f64 {
    let spin_variables = u32::from(spin_variables);
    let pow: f64 = 2f64.powi((spin_variables - 1) as i32);
    (PI.ln() * pow) - ln_gamma(pow)
}

/// Calculates the parameter complexity of an Independent Complete Component.
///
/// # Example
/// ```
/// # use mcm_finder_lib::mcm::parameter_complexity_icc;
/// # use std::num::NonZeroU32;
/// let spin_variables = NonZeroU32::new(9).unwrap();
/// assert_eq!(parameter_complexity_icc(spin_variables, 5), -58.3662038406751);
/// ```
pub fn parameter_complexity_icc(spin_variables: NonZeroU32, n: usize) -> f64 {
    let spin_variables = u32::from(spin_variables);
    let k = 2f64.powi(spin_variables as i32) - 1.0;
    k * (((n as f64) / 2f64) / PI).ln() / 2.0
}

/// The type for a Minimally Complex Model
///
/// # Example
/// ```
/// # use mcm_finder_lib::mcm::MinimallyComplexModel;
/// # use fixedbitset::FixedBitSet;
/// let mcm = MinimallyComplexModel::from_iccs(
///     vec![
///         FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
///         FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
///         FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
///     ],
/// ).unwrap();
/// assert_eq!(mcm.rank(), 9);
/// assert_eq!(mcm.complexity_mcm(), 1.3555732128424305);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MinimallyComplexModel {
    partition: Vec<IndependentCompleteComponent>,
}

#[derive(Debug, Clone, Copy)]
pub struct ParLogEResult {
    pub value: f64,
    pub new_icc_count: usize,
}

impl Deref for ParLogEResult {
    type Target = f64;
    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

impl MinimallyComplexModel {
    /// Create a new MCM with sorted ICCs.
    fn new_unsorted(partition: Vec<IndependentCompleteComponent>) -> MinimallyComplexModel {
        MinimallyComplexModel { partition }
    }

    /// Create a new MCM with sorted ICCs.
    fn new(partition: Vec<IndependentCompleteComponent>) -> MinimallyComplexModel {
        let mut mcm = MinimallyComplexModel::new_unsorted(partition);
        mcm.sort_first_bit();
        mcm
    }

    /// Returns an MCM from the given ICCs, with empty ICCs removed. Note that
    /// ICCs in an MCM will be sorted so may have a different order from how they were
    /// inserted.
    fn new_remove_empty(partition: Vec<IndependentCompleteComponent>) -> MinimallyComplexModel {
        let mut partition = partition;
        let variables = partition[0].len();
        partition.retain(|icc| icc.count_ones(..) > 0);
        if partition.is_empty() {
            partition.push(IndependentCompleteComponent::from(
                FixedBitSet::with_capacity_and_blocks(variables, [0]),
            ));
        }
        MinimallyComplexModel::new(partition)
    }

    /// Get the index of the ICC which contains this variable.
    fn get_index(&self, variable: usize) -> Option<usize> {
        for (nr, icc) in self.partition.iter().enumerate() {
            if icc[variable] {
                return Some(nr);
            }
        }
        None
    }

    /// Returns `true` if ICCs in array do not overlap anywhere. The ICCs do not
    /// have to include all of the variables.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mut partition = vec!(
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    /// );
    /// assert!(MinimallyComplexModel::verify_iccs(&partition));
    /// partition[1].set(0, true);
    /// assert!(!MinimallyComplexModel::verify_iccs(&partition));
    /// ```
    pub fn verify_iccs(partition: &[IndependentCompleteComponent]) -> bool {
        for (nr, one) in partition.iter().enumerate() {
            for two in &partition[(nr + 1)..] {
                if !(one & two).is_clear() {
                    return false;
                }
            }
        }
        true
    }

    /// Returns an empty MCM.
    pub fn empty(variables: NonZeroUsize) -> MinimallyComplexModel {
        MinimallyComplexModel {
            partition: vec![IndependentCompleteComponent::from(
                FixedBitSet::with_capacity_and_blocks(variables.into(), [0]),
            )],
        }
    }

    /// Sorts the ICCs in the MCM so that they look sorted by the first bit in
    /// the ICC when displayed.
    fn sort_first_bit(&mut self) {
        self.partition.sort_by(|a, b| {
            let a = a.as_slice().iter().map(|n| n.reverse_bits());
            let b = b.as_slice().iter().map(|n| n.reverse_bits());
            a.cmp(b).reverse()
        });
    }

    /// Create a new MCM from the given list of ICCs.
    ///
    /// This function returns `Ok(MinimallyComplexModel)` when the supplied ICCs
    /// have no overlap. Note that ICCs in an MCM will be sorted so may have a
    /// different order from how they were inserted.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mut partition = vec!(
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    /// );
    /// let mcm = MinimallyComplexModel::from_iccs(partition).unwrap();
    /// assert_eq!(mcm.rank(), 9);
    /// assert_eq!(mcm.count_icc(), 3);
    /// ```
    pub fn from_iccs(
        partition: Vec<icc::IndependentCompleteComponent>,
    ) -> Result<MinimallyComplexModel, MCMError> {
        if MinimallyComplexModel::verify_iccs(&partition) {
            return Ok(MinimallyComplexModel::new_remove_empty(partition));
        }
        Err(MCMError::FromIccs)
    }

    /// Create a new MCM from the given vector of numbers.
    ///
    /// The numbers in the vector should correspond to the ICC the variable of that
    /// index is a part of. '1' means that variable belongs to the first ICC. '0' means
    /// that variable is not included in the model.
    ///
    /// Note that ICCs in an MCM will be sorted so may have a different order from how
    /// they were inserted.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mut partition = vec!(
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    /// );
    /// let test_mcm = MinimallyComplexModel::from_iccs(partition).unwrap();
    /// assert_eq!(MinimallyComplexModel::from_vector(vec![1, 2, 2, 3, 3, 3, 2, 3, 3]), test_mcm);
    /// ```
    pub fn from_vector(vector: Vec<usize>) -> MinimallyComplexModel {
        let mut sorted = vector.clone();
        sorted.sort();
        sorted.dedup();
        sorted.retain(|&n| n > 0);
        let icc_amount = sorted.len();
        let map: Vec<(usize, &usize)> = sorted.iter().enumerate().collect();

        let continuous_vector: Vec<usize> = vector
            .into_iter()
            .map(|n| match n {
                0 => 0,
                1.. => map
                    .iter()
                    .find_map(|(nr, m)| if **m == n { Some(nr + 1) } else { None })
                    .unwrap(),
            })
            .collect();

        let mut iccs: Vec<icc::IndependentCompleteComponent> = vec![
            FixedBitSet::with_capacity_and_blocks(continuous_vector.len(), [0])
                .into();
            icc_amount
        ];
        for (nr, icc) in continuous_vector.into_iter().enumerate() {
            if icc > 0 {
                iccs[icc - 1].set(nr, true);
            }
        }

        MinimallyComplexModel::new(iccs)
    }

    /// Turn this MCM into a vector of integers representing the ICC each variable is
    /// included in. '1' means that variable belongs to the first ICC. '0' means that
    /// variable is not included in the model.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mut partition = vec!(
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    /// );
    /// let mcm = MinimallyComplexModel::from_iccs(partition).unwrap();
    /// assert_eq!(mcm.to_vector(), vec![1, 2, 2, 3, 3, 3, 2, 3, 3]);
    /// ```
    pub fn to_vector(&self) -> Vec<usize> {
        let mut output = vec![0usize; self.variables()];
        for (nr, icc) in self.partition.iter().enumerate() {
            for variable in icc.ones() {
                output[variable] = nr + 1;
            }
        }
        output
    }

    /// Returns an MCM with a single partition.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::full(NonZero::new(100).unwrap());
    /// assert_eq!(mcm.rank(), 100);
    /// assert_eq!(mcm.count_icc(), 1);
    /// ```
    pub fn full(variables: NonZeroUsize) -> MinimallyComplexModel {
        let mut partition = FixedBitSet::with_capacity(variables.into());
        partition.set_range(.., true);
        MinimallyComplexModel {
            partition: vec![partition.into()],
        }
    }

    /// Returns the trivial MCM with each variable in a different partition.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::trivial(NonZero::new(500).unwrap());
    /// assert_eq!(mcm.rank(), 500);
    /// assert_eq!(mcm.count_icc(), 500);
    /// ```
    pub fn trivial(variables: NonZeroUsize) -> MinimallyComplexModel {
        let variables: usize = variables.into();
        let mut partition = Vec::with_capacity(variables);
        let vector_length = variables.div_ceil(usize::BITS as usize);
        let mut part_content = vec![0usize; vector_length];
        let mut counter = 0usize;

        for n in 0..vector_length {
            part_content[n] = 1;
            for _ in 0..usize::BITS {
                partition.push(
                    FixedBitSet::with_capacity_and_blocks(variables, part_content.clone()).into(),
                );
                part_content[n] <<= 1;

                counter += 1;
                if counter == variables {
                    break;
                }
            }
        }
        MinimallyComplexModel { partition }
    }

    /// Returns the rank of this MCM. Can be lower than `self.variables()` if not
    /// all variables are included in the model.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    /// ]).unwrap();
    /// assert_eq!(mcm.rank(), 9);
    /// ```
    /// Now if we remove the first variable from this model:
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let smaller_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b010111000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    /// ]).unwrap();
    /// assert_eq!(smaller_mcm.rank(), 8);
    /// ```
    pub fn rank(&self) -> usize {
        self.partition
            .iter()
            .fold(FixedBitSet::with_capacity(self.variables()), |a, b| &a | b)
            .count_ones(..)
    }

    /// Returns the amount of variables in this MCM.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    /// ]).unwrap();
    /// assert_eq!(mcm.variables(), 9);
    /// ```
    pub fn variables(&self) -> usize {
        self.partition[0].len()
    }

    /// Merge the first basis into the second bases, and return a new MCM.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]).into(),
    /// ]).unwrap();
    /// assert_eq!(mcm.merge(2, 0), result_mcm);
    /// ```
    pub fn merge(&self, basis: usize, into: usize) -> MinimallyComplexModel {
        let mut iccs: Vec<icc::IndependentCompleteComponent> =
            self.partition.iter().map(|i| i.full_clone()).collect();
        // let mut iccs = self.partition.clone();

        iccs[into] |= &self.partition[basis];
        iccs.remove(basis);

        MinimallyComplexModel::new(iccs)
    }

    /// Split the marked variables from the basis ICC into a new ICC, and return a new MCM.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110101000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000010001]).into(),
    /// ]).unwrap();
    /// let mark = FixedBitSet::with_capacity_and_blocks(9, [0b001010001]);
    /// assert_eq!(mcm.split(0, mark), result_mcm);
    /// ```
    pub fn split(&self, basis: usize, split: FixedBitSet) -> MinimallyComplexModel {
        let mut mask = split.into();
        let mut iccs: Vec<icc::IndependentCompleteComponent> =
            self.partition.iter().map(|i| i.full_clone()).collect();
        // let mut iccs = self.partition.clone();
        let new_icc = &iccs[basis] & &mask;

        mask.toggle_range(..);
        iccs[basis] &= mask;

        iccs.push(new_icc);

        MinimallyComplexModel::new_remove_empty(iccs)
    }

    /// Split the marked variables from the basis ICC into a new ICC, and return a new MCM.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]).into(),
    /// ]).unwrap();
    /// println!("{}", mcm.split_one(0));
    /// assert_eq!(mcm.split_one(0), result_mcm);
    /// ```
    pub fn split_one(&self, variable: usize) -> MinimallyComplexModel {
        let basis = self.get_index(variable).unwrap();
        let mut iccs: Vec<icc::IndependentCompleteComponent> =
            self.partition.iter().map(|i| i.full_clone()).collect();
        // let mut iccs = self.partition.clone();
        let mut new_icc = FixedBitSet::with_capacity_and_blocks(iccs[0].len(), [0b0]);

        new_icc.set(variable, iccs[basis][variable]);
        iccs[basis].remove(variable);

        iccs.push(new_icc.into());

        MinimallyComplexModel::new_remove_empty(iccs)
    }

    /// Swaps the given variable from the basis ICC into the destination ICC.
    ///
    /// Swapping a variable into a non-existant ICC adds a new ICC.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]).into(),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001010110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110101001]).into(),
    /// ]).unwrap();
    /// assert_eq!(mcm.swap(4, 1), result_mcm);
    /// ```
    pub fn swap(&self, variable: usize, destination: usize) -> MinimallyComplexModel {
        let basis = self.get_index(variable);
        let mut iccs: Vec<IndependentCompleteComponent> =
            self.partition.iter().map(|i| i.full_clone()).collect();

        if destination >= self.count_icc() {
            return self.split_one(variable);
        }
        if let Some(basis) = basis {
            iccs[basis].set(variable, false);
        }
        iccs[destination].set(variable, true);

        MinimallyComplexModel::new_remove_empty(iccs)
    }

    /// Swaps the given variable from the basis ICC into the destination ICC.
    ///
    /// Swapping a variable into a non-existant ICC adds a new ICC. This method
    /// also returns the unchanged affected ICCs, for use in tabu search.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]).into(),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001010110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110101001]).into(),
    /// ]).unwrap();
    /// assert_eq!(mcm.swap_tabu(4, 1).0, result_mcm);
    /// ```
    pub fn swap_tabu(
        &self,
        variable: usize,
        destination: usize,
    ) -> (
        MinimallyComplexModel,
        Option<IndependentCompleteComponent>,
        Option<IndependentCompleteComponent>,
    ) {
        let basis = self.get_index(variable);
        if let Some(basis) = basis {
            // null operation
            if destination == basis {
                return (self.clone(), None, None);
            }

            // swap to new ICC
            if destination >= self.count_icc() {
                return (
                    self.split_one(variable),
                    Some(self.partition[basis].clone()),
                    None,
                );
            }

            let mut iccs: Vec<IndependentCompleteComponent> =
                self.partition.iter().map(|i| i.full_clone()).collect();

            iccs[basis].set(variable, false);
            iccs[destination].set(variable, true);
            (
                MinimallyComplexModel::new_remove_empty(iccs),
                Some(self.partition[basis].clone()),
                Some(self.partition[destination].clone()),
            )
        } else {
            // missing variable
            let mut iccs: Vec<IndependentCompleteComponent> =
                self.partition.iter().map(|i| i.full_clone()).collect();
            iccs[destination].set(variable, true);

            (
                MinimallyComplexModel::new_remove_empty(iccs),
                None,
                Some(self.partition[destination].clone()),
            )
        }
    }

    /// Combine the two ICCs, then put the unmasked variables in one ICC and the
    /// masked variables in the other.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]).into(),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001010110]).into(),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110101001]).into(),
    /// ]).unwrap();
    /// let mask = dbg!(FixedBitSet::with_capacity_and_blocks(9, [0b001010110]));
    /// let result = dbg!(mcm.redistribute(0, 1, mask));
    /// assert_eq!(result, result_mcm);
    /// ```
    pub fn redistribute(
        &self,
        icc_one: usize,
        icc_two: usize,
        mask: FixedBitSet,
    ) -> MinimallyComplexModel {
        assert_ne!(icc_one, icc_two);
        let new_mcm = self.merge(icc_two, icc_one);
        new_mcm.split(icc_one.min(icc_two), mask)
    }

    /// Returns the amount of ICCs present in this model.
    ///
    /// # Example
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::trivial(NonZero::new(5).unwrap());
    /// assert_eq!(mcm.count_icc(), 5);
    /// ```
    pub fn count_icc(&self) -> usize {
        self.partition.len()
    }

    pub fn is_tabu(&self, tabu_list: &VecDeque<IndependentCompleteComponent>) -> bool {
        self.partition.iter().any(|icc| tabu_list.contains(icc))
    }

    pub fn complexity_mcm(&self) -> f64 {
        let mut c_param = 0.0f64;
        let mut c_geom = 0.0f64;

        for part in self.partition.iter() {
            let spin_variables = NonZeroU32::try_from(part.count_ones(..) as u32).unwrap();
            c_param += parameter_complexity_icc(spin_variables, self.rank());
            c_geom += geometric_complexity_icc(spin_variables);
        }

        c_param + c_geom
    }

    /// Returns the amount of ICCs in the model with more than one variable.
    ///
    /// # Example
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::trivial(NonZero::new(5).unwrap());
    /// assert_eq!(mcm.merge(2, 3).count_nontrivial_icc(), 1);
    /// ```
    pub fn count_nontrivial_icc(&self) -> usize {
        self.partition
            .iter()
            .filter(|icc| icc.count_ones(..) > 1)
            .count()
    }

    /// Returns a new MCM with a random ICC mutation
    pub fn mutate(&self, rng: &mut rand::rngs::ThreadRng) -> (MinimallyComplexModel, MutationType) {
        let mut_type = if self.count_icc() == 1 {
            MutationType::Split
        } else if self.count_nontrivial_icc() == 0 {
            MutationType::Merge
        } else {
            MutationType::rand(rng, [1.0, 1.0, 1.0])
        };

        let new_mcm = match mut_type {
            MutationType::Merge => {
                let targets = (0..self.count_icc()).sample(rng, 2);
                self.merge(targets[0], targets[1])
            }
            MutationType::Split => {
                let candidates: Vec<usize> = self
                    .partition
                    .iter()
                    .enumerate()
                    .filter_map(|(nr, icc)| {
                        if icc.count_ones(..) > 1 {
                            Some(nr)
                        } else {
                            None
                        }
                    })
                    .collect();
                let basis = *candidates.choose(rng).unwrap();
                // let mut random_data: Vec<u64> = vec![0; self.variables().div_ceil(32)];
                // rng.fill(&mut random_data);

                let variable_count = self.partition[basis].count_ones(..);
                let distribution = Binomial::new(0.5, (variable_count - 2) as u64).unwrap();
                let sample: u64 = distribution.sample(rng);
                // let amount = rng.random_range(1..variable_count);
                let amount = (sample + 1) as usize;
                let random_variables = self.partition[basis].ones().sample(rng, amount);

                let mut split = FixedBitSet::with_capacity(self.variables());
                for i in random_variables {
                    split.insert(i);
                }

                // let split = FixedBitSet::with_capacity_and_blocks(
                //     self.variables(),
                //     random_data.into_iter().map(|n| n as usize),
                // );
                self.split(basis, split)
                // let choice = self.partition[basis].ones().choose(rng).unwrap();
                // self.split_one(basis, choice)
            }
            MutationType::Swap => {
                // let basis_candidates: Vec<usize> = self
                //     .partition
                //     .iter()
                //     .enumerate()
                //     .filter_map(|(nr, icc)| {
                //         if icc.count_ones(..) > 1 {
                //             Some(nr)
                //         } else {
                //             None
                //         }
                //     })
                //     .collect();
                // let basis = *basis_candidates.choose(rng).unwrap();
                // let destination = (0usize..self.count_icc())
                //     .filter(|n| *n != basis)
                //     .choose(rng)
                //     .unwrap();
                // let choice = self.partition[basis].ones().choose(rng).unwrap();

                let (choice, candidate_index) = loop {
                    // This will loop forever when the MCM is empty, but don't mutate an empty MCM.
                    let candidate_variable = (0..self.variables()).choose(rng).unwrap_or_default();
                    if let Some(index) = self.get_index(candidate_variable) {
                        break (candidate_variable, index);
                    }
                };
                let destination = (0usize..=self.count_icc())
                    .filter(|n| *n != candidate_index)
                    .choose(rng)
                    .unwrap();

                self.swap(choice, destination)
            }
        };
        (new_mcm, mut_type)
    }

    /// Calculate the logarithm of the evidence of this MCM, via the equation
    ///
    /// log⁡E=-(n-r)ln2 + ∑_{P∈P}\[(ln⁡Γ(2ᴾ⁻¹)−ln⁡Γ(n+2ᴾ⁻¹))+∑_{x∈data}(P)(ln⁡Γ(kx+0.5)−ln⁡Γ(0.5))\]
    ///
    /// Hope that clears things up.
    pub fn log_e<T: Dataset>(
        &self,
        dataset: &T,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
    ) -> f64 {
        let mut log_e = 0f64;

        for icc in self.partition.iter() {
            if let Some(cache) = log_e_cache {
                if let Some(icc_log_e) = icc.log_e.get() {
                    log_e += icc_log_e;
                } else {
                    log_e += *cache
                        .entry(icc.bits.clone())
                        .or_insert_with(|| icc.log_e(dataset));
                }
            } else {
                log_e += icc.log_e(dataset);
            }
        }

        let front_constant: f64 =
            (dataset.datapoints() * (self.variables() - self.rank())) as f64 * LN_2;
        log_e - front_constant
    }

    /// Calculate the logarithm of the evidence of this MCM, parallelized with Rayon.
    pub fn par_log_e<T: Dataset + Sync>(
        &self,
        dataset: &T,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
    ) -> ParLogEResult {
        let log_e: (f64, usize) = self
            .partition
            .par_iter()
            .map(|icc| {
                if let Some(cache) = log_e_cache {
                    if let Some(icc_log_e) = icc.log_e.get() {
                        (*icc_log_e, 0)
                    } else {
                        let mut hit = 1;
                        let result = *cache.entry(icc.bits.clone()).or_insert_with(|| {
                            hit -= 1;
                            icc.log_e(dataset)
                        });
                        (result, hit)
                    }
                } else {
                    (icc.log_e(dataset), 1)
                }
            })
            .reduce(|| (0.0, 0), |a, b| (a.0 + b.0, a.1 + b.1));

        let front_constant: f64 =
            (dataset.datapoints() * (self.variables() - self.rank())) as f64 * LN_2;
        ParLogEResult {
            value: log_e.0 - front_constant,
            new_icc_count: log_e.1,
        }
    }

    pub fn read_from_file(path: &Path) -> Result<MinimallyComplexModel, MCMError> {
        let mut iccs: Vec<icc::IndependentCompleteComponent> = vec![];
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

                    // add the datapoints to the vector
                    let mut bitvec = FixedBitSet::with_capacity(bool_array.len());
                    for (nr, bit) in bool_array.drain(..).enumerate() {
                        bitvec.set(nr, bit);
                    }
                    iccs.push(bitvec.into());
                    debug_assert!(bool_array.is_empty());
                }
                b'\r' | b'\n' => {}
                // wrong character case
                _ => Err(MCMError::BadCharacter {
                    src: NamedSource::new(&filename, file.clone()),
                    bad_line: nr.into(),
                })?,
            }
        }

        Ok(MinimallyComplexModel { partition: iccs })
    }

    pub fn to_matrix(&self) -> Vec<FixedBitSet> {
        let mut output =
            vec![FixedBitSet::with_capacity_and_blocks(self.variables(), [0]); self.variables()];

        for icc in self.partition.iter() {
            for row in icc.bits.ones() {
                for col in icc.bits.ones() {
                    output[row].set(col, true);
                }
            }
        }

        output
    }

    pub fn jakkard_difference(&self, other: &MinimallyComplexModel) -> f64 {
        let self_matrix = self.to_matrix();
        let other_matrix = other.to_matrix();

        let mut or_total = 0usize;
        let mut and_total = 0usize;
        for (s, o) in self_matrix.iter().zip(other_matrix.iter()) {
            or_total += (s | o).count_ones(..);
            and_total += (s & o).count_ones(..);
        }

        and_total as f64 / or_total as f64
    }
}

fn gamma_factor<T: Dataset>(dataset: &T, rank_subset: i32) -> f64 {
    let points = dataset.datapoints() as f64;
    if rank_subset > 25 {
        (-points * ((rank_subset - 1) as f64) * LN_2)
            - (points * (points - 1.0)) / (2.0f64.powi(rank_subset))
    } else {
        ln_gamma(2.0f64.powi(rank_subset - 1)) - ln_gamma(points + 2.0f64.powi(rank_subset - 1))
    }
}

impl Display for MinimallyComplexModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut partition_string: String = self
            .partition
            .iter()
            .flat_map(|p| format!("{p}\n").chars().collect::<Vec<_>>())
            .collect();
        partition_string.pop();
        write!(f, "{}", partition_string)
    }
}
