use std::{collections::HashMap, fmt::Display, marker::Sized, path::Path};

use crate::{mcm::MinimallyComplexModel, mcm_error::MCMError};

/// A class contaning the result of a `Solver`. Contains:
///  - `mcm`: The MCM found by this algorithm with the highest Log(E)
///  - `log_e`: The Log(E) of this MCM
///  - `other_stuff`: other information relevant to this result
#[derive(Debug, Clone)]
pub struct SolverReport {
    pub mcm: MinimallyComplexModel,
    pub log_e: f64,
    pub other_stuff: HashMap<String, String>,
}

impl SolverReport {
    /// Create a new `SolverReport`.
    ///
    /// `other_stuff` can be used to return miscellaneous data about the result
    /// or the sorting process.
    pub(crate) fn new(
        mcm: MinimallyComplexModel,
        log_e: f64,
        other_stuff: HashMap<String, String>,
    ) -> SolverReport {
        SolverReport {
            mcm,
            log_e,
            other_stuff,
        }
    }
}

impl Display for SolverReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Found MCM with Log E of {} and ICCs:\n{}{}",
            self.log_e,
            self.mcm,
            self.other_stuff
                .iter()
                .flat_map(|(k, v)| format!("\n - {:?}: {:?}", k, v).chars().collect::<Vec<_>>())
                .collect::<String>()
        )
    }
}

pub trait Solver {
    /// Load data from the given path.
    ///
    /// Returns a `Result` containing the `Solver` or an `MCMError` if the data
    /// reading or loading went wrong.
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized;

    /// Start the solver.
    ///
    /// Returns a `SolverReport` which can be printed or used further.
    fn solve(&self) -> SolverReport;
}
