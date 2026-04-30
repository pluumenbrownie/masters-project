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
}

impl Solver for GreedySearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(GreedySearcher {
            dataset: VecDataset::read_from_file(filepath)?,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut best_mcm =
            MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap());

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
            let mut new_best = false;
            let old_best = best_mcm.clone();
            progress.set_description(format!("{parts_left} steps left"));
            for basis in 1..=parts_left {
                for into in 0..basis {
                    let candidate = old_best.merge(basis, into);
                    if candidate
                        .log_e(&self.dataset, &mut log_e_cache)
                        .total_cmp(&best_log_e)
                        .is_ge()
                    {
                        best_mcm = candidate.clone();
                        best_log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache);
                        new_best = true;
                    }
                    let _ = progress.update(1);
                }
            }
            if !new_best {
                break;
            }
        }

        SolverReport::new(best_mcm, best_log_e, HashMap::new())
    }
}
