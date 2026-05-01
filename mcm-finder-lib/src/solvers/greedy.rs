use std::{collections::HashMap, marker::Sized, num::NonZero, path::Path};

use kdam::{BarExt, tqdm};

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::solvers_base::{Solver, SolverReport},
};

pub struct GreedySearcher {
    pub(crate) dataset: VecDataset,
    continue_after_minimum: bool,
}

impl GreedySearcher {
    /// Toggle whether the greedy algorithm should stop searching when a step does not
    /// generate a new minimum (default) or should finish all steps.
    pub fn continue_after_minimum(mut self) -> Self {
        self.continue_after_minimum = true;
        self
    }
}

impl Solver for GreedySearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(GreedySearcher {
            dataset: VecDataset::read_from_file(filepath)?,
            continue_after_minimum: false,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut best_mcm =
            MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap());
        let mut gen_best = best_mcm.clone();

        // let mut best_mcm: MinimallyComplexModel = current.clone();
        let mut log_e_cache = Some(HashMap::new());
        let mut best_log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache);

        let mut length = 0usize;
        for parts_left in (0usize..self.dataset.variables()).rev() {
            for basis in tqdm!(1..=parts_left) {
                for _ in 0..basis {
                    length += 1;
                }
            }
        }

        // we merge one partition each round
        let mut progress = tqdm!(total = length);
        for parts_left in (0usize..self.dataset.variables()).rev() {
            let old_best = gen_best.clone();
            let mut new_best = false;
            progress.set_description(format!("{parts_left} steps left"));
            for basis in 1..=parts_left {
                for into in 0..basis {
                    let candidate = old_best.merge(basis, into);

                    // update current best
                    if candidate
                        .log_e(&self.dataset, &mut log_e_cache)
                        .total_cmp(&gen_best.log_e(&self.dataset, &mut log_e_cache))
                        .is_ge()
                    {
                        gen_best = candidate.clone();
                    }

                    // update global best
                    if candidate
                        .log_e(&self.dataset, &mut log_e_cache)
                        .total_cmp(&best_mcm.log_e(&self.dataset, &mut log_e_cache))
                        .is_ge()
                    {
                        best_mcm = candidate.clone();
                        best_log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache);
                        new_best = true;
                    }
                    let _ = progress.update(1);
                }
            }
            if !new_best & !self.continue_after_minimum {
                break;
            }
        }

        SolverReport::new(best_mcm, best_log_e, HashMap::new())
    }
}
