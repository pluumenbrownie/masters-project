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

    pub fn set_deep_log_e(self, deep_log_e: f64) -> Self {
        Self {
            mcm: self.mcm,
            log_e: self.log_e,
            deep_log_e,
        }
    }
}

#[derive(Clone)]
pub struct GreedySearcher {
    dataset: VecDataset,
    lookahead_depth: usize,
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
    ///
    /// Setting this value will explode the execution time of this algorithm, so use
    /// with caution.
    pub fn lookahead(mut self, depth: usize) -> Self {
        self.lookahead_depth = depth;
        self
    }

    fn update_global_best(&self, best_mcm: &mut LogeMCM, candidate: &LogeMCM) -> bool {
        if candidate.log_e.total_cmp(&best_mcm.log_e).is_gt() {
            *best_mcm = candidate.clone();
            return true;
        };
        false
    }

    #[must_use]
    fn find_best_merge(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        progress: &mut Bar,
        iccs_left: usize,
        original: &LogeMCM,
        depth: usize,
    ) -> Option<LogeMCM> {
        if depth > self.lookahead_depth || iccs_left == 0 {
            None
        } else {
            let mut gen_best: Option<LogeMCM> = None;
            for basis in 1..=iccs_left {
                for into in 0..basis {
                    let mut candidate = LogeMCM::calculate(
                        original.mcm.merge(basis, into),
                        &self.dataset,
                        log_e_cache,
                    );

                    // lookahead
                    if let Some(deep_best) = self.find_best_merge(
                        log_e_cache,
                        progress,
                        iccs_left - 1,
                        &candidate,
                        depth + 1,
                    ) {
                        let deep_log_e = deep_best.deep_log_e;
                        candidate = update_deep_log_e_only(candidate, deep_log_e);
                    }

                    gen_best = update_if_deep_log_e_better(gen_best, &candidate);

                    let _ = progress.update(1);
                }
            }
            gen_best
        }
    }

    fn count_calculations(&self) -> usize {
        let mut length = 0usize;
        let depth = self.lookahead_depth;

        for parts_left in (0usize..self.dataset.variables()).rev() {
            count_calculations_recursive(&mut length, parts_left, depth);
        }
        length
    }
}
fn count_calculations_recursive(length: &mut usize, parts_left: usize, depth: usize) {
    for basis in 1..=parts_left {
        for _ in 0..basis {
            *length += 1;
            if depth > 0 {
                count_calculations_recursive(length, parts_left - 1, depth - 1);
            }
        }
    }
}

fn update_deep_log_e_only(candidate: LogeMCM, deep_log_e: f64) -> LogeMCM {
    if deep_log_e.total_cmp(&candidate.deep_log_e).is_gt() {
        candidate.set_deep_log_e(deep_log_e)
    } else {
        candidate
    }
}

fn update_if_deep_log_e_better(gen_best: Option<LogeMCM>, candidate: &LogeMCM) -> Option<LogeMCM> {
    gen_best
        .filter(|gen_best| gen_best.deep_log_e.total_cmp(&candidate.deep_log_e).is_ge())
        .or(Some(candidate.clone()))
}

impl Solver for GreedySearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(GreedySearcher {
            dataset: VecDataset::read_from_file(filepath)?,
            lookahead_depth: 0,
            continue_after_minimum: false,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut log_e_cache = get_log_e_cache();

        let mut best_mcm = LogeMCM::calculate(
            MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap()),
            &self.dataset,
            &mut log_e_cache,
        );
        let mut gen_best = best_mcm.clone();

        let length = self.count_calculations();

        // we merge one partition each round
        let mut progress = tqdm!(total = length);
        for iccs_left in (0usize..self.dataset.variables()).rev() {
            let original = gen_best.clone();
            progress.set_description(format!("{iccs_left} ICCs"));
            gen_best = self
                .find_best_merge(&mut log_e_cache, &mut progress, iccs_left, &original, 0)
                .unwrap();
            let new_best = self.update_global_best(&mut best_mcm, &gen_best);
            if !new_best & !self.continue_after_minimum {
                break;
            }
        }

        SolverReport::new(
            best_mcm.mcm,
            best_mcm.log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}
