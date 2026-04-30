use std::{collections::HashMap, fmt::Display, marker::Sized, path::Path};

use crate::{mcm::MinimallyComplexModel, mcm_error::MCMError};

#[derive(Debug, Clone)]
pub struct SolverReport {
    pub mcm: MinimallyComplexModel,
    pub log_e: f64,
    pub other_stuff: HashMap<String, String>,
}

impl SolverReport {
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
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized;

    fn solve(&self) -> SolverReport;
}
