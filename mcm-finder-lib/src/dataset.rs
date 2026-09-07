//! This module contains the different kinds of datasets and relevant traits in
//! this library.

use std::{hash::RandomState, path::Path};

use fixedbitset::FixedBitSet;
use itertools::Itertools;
use miette::NamedSource;

use crate::mcm_error::MCMError;
pub use ends_cached::*;
pub use simple::*;

pub mod datacontainer;
mod ends;
pub mod ends_cached;
pub mod simple;

/// The top level trait for datasets. Datasets can:
///
/// - Read data from a text file.
/// - Return basic properties of the loaded data; and
/// - Calculate the log evidence of the data for a given ICC;
///
/// The [`Dataset`] trait is agnostic to the underlying way the data is stored,
/// or how the log evidence is calculated.
pub trait Dataset {
    /// Returns the amount of datapoints in this dataset.
    fn datapoints(&self) -> usize;

    /// Returns the amount of variables in each datapoint.
    fn variables(&self) -> usize;

    /// Returns the amount of states each variable can have. This count is detected
    /// automatically in the `read_from_file` method.
    fn variable_states(&self) -> usize;

    /// Returns the bit width of a variable.
    fn variable_width(&self) -> usize {
        (self.variable_states() as f64).log2().ceil() as usize
    }

    /// Returns the prevalence of each state of this variable.
    fn state_prevalence(&self, variable: usize) -> Vec<usize>;

    /// Enlarge the variables in this ICC to cover multi-bit variables.
    fn resize_icc(&self, icc: &FixedBitSet) -> FixedBitSet {
        resize_mask(icc.clone(), self.variable_width())
    }

    /// Compute the logarithmic evidence for the given ICC.
    fn log_e(&self, icc: &FixedBitSet) -> f64;

    fn read_from_file(path: &Path) -> Result<Self, MCMError>
    where
        Self: Sized;
}

pub type DefaultState = std::hash::RandomState;
pub type AhashState = ahash::RandomState;
pub type RapidStateFast = rapidhash::fast::RandomState;
pub type RapidStateQuality = rapidhash::quality::RandomState;
pub type FxState = rustc_hash::FxBuildHasher;

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
            help_line: Some(format!(
                "Line has length of {}, instead of detected line length {}.",
                bool_array.len(),
                line_length
            )),
        })?
    };
    Ok(())
}

/// Enlarge the variables in this bitmask to cover multi-bit variables.
///
/// # Examples
/// ```
/// # use fixedbitset::FixedBitSet;
/// # use mcm_finder_lib::dataset::resize_mask;
///
/// let small_mask = FixedBitSet::with_capacity_and_blocks(5, [0b01001]);
/// let big_mask = FixedBitSet::with_capacity_and_blocks(15, [0b000111000000111]);
///
/// let result = resize_mask(small_mask, 3);
///
/// assert_eq!(result, big_mask);
/// ```
pub fn resize_mask(mask: FixedBitSet, variable_width: usize) -> FixedBitSet {
    if variable_width == 1 {
        return mask;
    }
    let mut big_mask = FixedBitSet::with_capacity(mask.len() * variable_width);
    for bit in mask.ones() {
        for i in 0..variable_width {
            big_mask.set(bit * variable_width + i, true);
        }
    }
    big_mask
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

pub(crate) fn get_and_check_unique_symbols(
    filename: &str,
    file: &str,
) -> Result<Vec<u8>, MCMError> {
    let unique_symbols = file.bytes().enumerate().unique_by(|s| s.1).collect_vec();
    unique_symbols
        .iter()
        .map(|(nr, b)| (nr, b.is_ascii_alphanumeric() || *b == b'\r' || *b == b'\n'))
        .map(|(&nr, b)| {
            b.ok_or_else(|| MCMError::BadCharacter {
                src: NamedSource::new(filename, file.to_string()),
                bad_line: nr.into(),
            })
        })
        .find(|r| r.is_err())
        .unwrap_or(Ok(()))?;
    Ok(unique_symbols
        .into_iter()
        .map(|(_, b)| b)
        .filter(|b| b.is_ascii_alphanumeric())
        .sorted()
        .collect_vec())
}
