use std::{
    collections::HashMap,
    marker::Sized,
    num::{NonZero, NonZeroUsize},
    path::Path,
};

use fixedbitset::FixedBitSet;
use kdam::{Bar, BarExt, tqdm};

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::solvers_base::{Solver, SolverReport},
};

#[derive(Debug, Clone)]
pub struct LogeMCM {
    pub mcm: MinimallyComplexModel,
    pub log_e: f64,
    pub deep_log_e: f64,
}

impl LogeMCM {
    pub fn new(mcm: MinimallyComplexModel, log_e: f64) -> Self {
        Self {
            mcm,
            log_e,
            deep_log_e: log_e,
        }
    }

    pub fn calculate<T: Dataset>(
        mcm: MinimallyComplexModel,
        dataset: &T,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
    ) -> Self {
        let log_e = mcm.log_e(dataset, log_e_cache);
        Self {
            mcm,
            log_e,
            deep_log_e: log_e,
        }
    }
}

#[derive(Clone)]
pub struct GreedySearcher {
    dataset: VecDataset,
    lookahead: Option<NonZeroUsize>,
    continue_after_minimum: bool,
}

impl GreedySearcher {
    /// Toggle whether the greedy algorithm should stop searching when a step does not
    /// generate a new minimum (default) or should finish all steps.
    pub fn continue_after_minimum(mut self) -> Self {
        self.continue_after_minimum = true;
        self
    }

    /// Sets the amount of steps the algorithm is allowed to look ahead to choose the
    /// next best candidate.
    pub fn lookahead(mut self, depth: usize) -> Self {
        self.lookahead =
            Some(NonZeroUsize::try_from(depth).expect("Lookahead depth should not be zero."));
        self
    }

    fn update_global_best(&self, best_mcm: &mut LogeMCM, candidate: &LogeMCM) -> bool {
        if candidate.log_e.total_cmp(&best_mcm.log_e).is_ge() {
            *best_mcm = candidate.clone();
            return true;
        };
        false
    }

    #[must_use]
    fn find_lowest_merge(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        progress: &mut Bar,
        parts_left: usize,
        old_best: &LogeMCM,
    ) -> LogeMCM {
        let mut gen_best: Option<LogeMCM> = None;
        for basis in 1..=parts_left {
            for into in 0..basis {
                let candidate =
                    LogeMCM::calculate(old_best.mcm.merge(basis, into), &self.dataset, log_e_cache);

                match &gen_best {
                    Some(best) => {
                        if candidate.deep_log_e.total_cmp(&best.deep_log_e).is_ge() {
                            gen_best = Some(candidate.clone());
                        }
                    }
                    None => gen_best = Some(candidate.clone()),
                }

                let _ = progress.update(1);
            }
        }
        gen_best.unwrap()
    }

    fn count_calculations(&self) -> usize {
        let mut length = 0usize;
        for parts_left in (0usize..self.dataset.variables()).rev() {
            for basis in 1..=parts_left {
                for _ in 0..basis {
                    length += 1;
                }
            }
        }
        length
    }
}

impl Solver for GreedySearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(GreedySearcher {
            dataset: VecDataset::read_from_file(filepath)?,
            lookahead: None,
            continue_after_minimum: false,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut log_e_cache = Some(HashMap::new());

        let mut best_mcm = LogeMCM::calculate(
            MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap()),
            &self.dataset,
            &mut log_e_cache,
        );
        let mut gen_best = best_mcm.clone();

        let length = self.count_calculations();

        // we merge one partition each round
        let mut progress = tqdm!(total = length);
        for parts_left in (0usize..self.dataset.variables()).rev() {
            let old_best = gen_best.clone();
            progress.set_description(format!("{parts_left} steps left"));
            gen_best =
                self.find_lowest_merge(&mut log_e_cache, &mut progress, parts_left, &old_best);
            let new_best = self.update_global_best(&mut best_mcm, &gen_best);
            if !new_best & !self.continue_after_minimum {
                break;
            }
        }

        SolverReport::new(best_mcm.mcm, best_mcm.log_e, HashMap::new())
    }
}
