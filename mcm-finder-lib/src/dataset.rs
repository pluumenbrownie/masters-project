//! This module contains the different kinds of datasets and relevant traits in
//! this library.

use std::{hash::RandomState, path::Path};

use fixedbitset::FixedBitSet;
use miette::NamedSource;

use crate::mcm_error::MCMError;

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

    /// Returns the prevalence of each state of this variable.
    fn state_prevalence(&self, variable: usize) -> Vec<usize>;

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
