use std::{collections::HashMap, marker::Sized, num::NonZero, path::Path};

use fixedbitset::FixedBitSet;
use kdam::{Bar, BarExt, tqdm};

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::{
        get_log_e_cache,
        solvers_base::{Solver, SolverReport},
    },
};

#[derive(Clone)]
pub struct ConstructiveSolver {
    dataset: VecDataset,
}

impl Solver for ConstructiveSolver {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(ConstructiveSolver {
            dataset: VecDataset::read_from_file(filepath)?,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut log_e_cache = get_log_e_cache();

        let mut current_solution =
            MinimallyComplexModel::empty(self.dataset.variables().try_into().unwrap());

        for var in tqdm!(0..self.dataset.variables()) {
            let vec_rep = current_solution.to_vector();
            let max_icc = vec_rep.iter().max().unwrap() + 1;

            let candidates: Vec<(MinimallyComplexModel, f64)> = (1..=max_icc)
                .map(|c| {
                    let mut new_candidate = vec_rep.clone();
                    new_candidate[var] = c;
                    let mcm = MinimallyComplexModel::from_vector(new_candidate);
                    let log_e = mcm.log_e(&self.dataset, &mut log_e_cache);
                    (mcm, log_e)
                })
                .collect();

            current_solution = candidates
                .into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap()
                .0;
        }

        let log_e = current_solution.log_e(&self.dataset, &mut log_e_cache);
        SolverReport::new(
            current_solution,
            log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}
