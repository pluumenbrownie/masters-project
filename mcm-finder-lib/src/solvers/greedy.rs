use std::{
    char::ToLowercase, collections::HashMap, default, marker::Sized, num::NonZero,
    ops::RangeBounds, path::Path,
};

use fixedbitset::FixedBitSet;
use kdam::{Bar, BarExt, tqdm};
use rand::seq::{IteratorRandom, SliceRandom};

use crate::{
    dataset::{Dataset, simple::VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::{
        ConstructiveStrategy, get_log_e_cache,
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

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum InitialSolver {
    #[default]
    Merge,
    MergeBeam {
        amount: usize,
    },
    Construct,
    None,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Refinements {
    Local,
    ChooseN { n: usize, max_fails: usize },
}

#[derive(Clone)]
pub struct GreedySolver {
    dataset: VecDataset,
    continue_after_minimum: bool,
    initial_solver: InitialSolver,
    constructive_strategy: ConstructiveStrategy,
    refinement_sequence: Vec<Refinements>,
}

impl GreedySolver {
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
    pub fn set_initial_solver(mut self, solver: InitialSolver) -> Self {
        self.initial_solver = solver;
        self
    }

    /// What search algorithms should be used after the initial construction algorithm.
    pub fn set_refinement_sequence(mut self, sequence: Vec<Refinements>) -> Self {
        self.refinement_sequence = sequence;
        self
    }

    fn update_global_best(&self, best_mcm: &mut LogeMCM, candidate: &LogeMCM) -> bool {
        if candidate.log_e.total_cmp(&best_mcm.log_e).is_gt() {
            *best_mcm = candidate.clone();
            return true;
        };
        false
    }

    // #[must_use]
    // fn find_best_merge(
    //     &self,
    //     log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
    //     progress: &mut Bar,
    //     iccs_left: usize,
    //     original: &LogeMCM,
    //     depth: usize,
    // ) -> Option<LogeMCM> {
    //     if depth > self.lookahead_depth || iccs_left == 0 {
    //         None
    //     } else {
    //         let mut gen_best: Option<LogeMCM> = None;
    //         for basis in 1..=iccs_left {
    //             for into in 0..basis {
    //                 let mut candidate = LogeMCM::calculate(
    //                     original.mcm.merge(basis, into),
    //                     &self.dataset,
    //                     log_e_cache,
    //                 );

    //                 // lookahead
    //                 if let Some(deep_best) = self.find_best_merge(
    //                     log_e_cache,
    //                     progress,
    //                     iccs_left - 1,
    //                     &candidate,
    //                     depth + 1,
    //                 ) {
    //                     let deep_log_e = deep_best.deep_log_e;
    //                     candidate = update_deep_log_e_only(candidate, deep_log_e);
    //                 }

    //                 gen_best = update_if_deep_log_e_better(gen_best, &candidate);

    //                 let _ = progress.update(1);
    //             }
    //         }
    //         gen_best
    //     }
    // }

    #[must_use]
    fn find_best_merge(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        progress: &mut Bar,
        iccs_left: usize,
        original_vec: &[LogeMCM],
        beam_size: usize,
    ) -> Vec<LogeMCM> {
        if iccs_left == 0 {
            return Vec::new();
        }

        let mut gen_best: Vec<LogeMCM> = vec![];
        for original in original_vec {
            for basis in 1..=iccs_left {
                for into in 0..basis {
                    let candidate = LogeMCM::calculate(
                        original.mcm.merge(basis, into),
                        &self.dataset,
                        log_e_cache,
                    );

                    if gen_best
                        .last()
                        .is_none_or(|o| o.log_e.total_cmp(&candidate.log_e).is_lt())
                        & gen_best.iter().all(|o| o.mcm != candidate.mcm)
                    {
                        if gen_best.len() >= beam_size {
                            gen_best.pop();
                        }
                        gen_best.push(candidate);
                        gen_best.sort_by(|a, b| a.log_e.total_cmp(&b.log_e).reverse());
                    }

                    let _ = progress.update(1);
                }
            }
        }
        gen_best
    }

    #[must_use]
    fn find_best_local(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        progress: &mut Bar,
        original: &LogeMCM,
    ) -> Option<LogeMCM> {
        let icc_amount = original.mcm.count_icc();

        let mut gen_best = None;
        for variable in 0..original.mcm.variables() {
            for into in 0..=icc_amount {
                let candidate = LogeMCM::calculate(
                    original.mcm.swap(variable, into),
                    &self.dataset,
                    log_e_cache,
                );

                gen_best = update_if_deep_log_e_better(gen_best, &candidate);
                let _ = progress.update(1);
            }
        }

        gen_best
    }

    #[must_use]
    fn find_best_check_amount(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        progress: &mut Bar,
        amount: usize,
        original: &LogeMCM,
    ) -> Option<LogeMCM> {
        let mut rng = rand::rng();
        let icc_amount = original.mcm.count_icc();

        let selection = (0..original.mcm.variables()).sample(&mut rng, amount);
        let icc_iterator = ProductIterator::new(amount, icc_amount);
        let total_length = icc_iterator.len();
        progress.total += total_length;

        let mut gen_best = None;
        for into_vec in icc_iterator {
            let mut mcm = original.mcm.clone();
            for (var, into) in selection.iter().zip(into_vec.iter()) {
                mcm = mcm.swap(*var, *into);
            }
            let candidate = LogeMCM::calculate(mcm, &self.dataset, log_e_cache);

            gen_best = update_if_deep_log_e_better(gen_best, &candidate);
            let _ = progress.update(1);
        }

        gen_best
    }

    fn count_calculations(&self) -> usize {
        let mut length = 0usize;

        for parts_left in (0usize..self.dataset.variables()).rev() {
            count_calculations(&mut length, parts_left);
        }
        length
    }

    fn solve_merge(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        best_mcm: &mut LogeMCM,
        amount: usize,
    ) -> Bar {
        let mut gen_best_vec = vec![best_mcm.clone()];
        let length = self.count_calculations();

        // we merge one partition each round
        let mut progress = tqdm!(total = length * amount);
        for iccs_left in (1usize..self.dataset.variables()).rev() {
            let original_vec = gen_best_vec.clone();
            progress.set_description(format!("{iccs_left} ICCs - {:.0}", gen_best_vec[0].log_e));
            gen_best_vec =
                self.find_best_merge(log_e_cache, &mut progress, iccs_left, &original_vec, amount);
            let new_best = self.update_global_best(best_mcm, &gen_best_vec[0]);
            if !new_best & !self.continue_after_minimum {
                break;
            }
            progress.update(1).unwrap();
        }
        progress
    }

    fn solve_local(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        best_mcm: &mut LogeMCM,
        progress: &mut Bar,
    ) {
        if self.initial_solver != InitialSolver::None {
            println!("{}", best_mcm.mcm);
        }
        loop {
            let original = best_mcm.clone();
            progress.set_description(format!("Local - {:.0}", best_mcm.log_e));
            if self
                .find_best_local(log_e_cache, progress, &original)
                .and_then(|c| self.update_global_best(best_mcm, &c).then_some(()))
                .is_none()
            {
                break;
            };
        }
    }

    fn solve_check_amount(
        &self,
        amount: usize,
        max_fails: usize,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        best_mcm: &mut LogeMCM,
        progress: &mut Bar,
    ) {
        let mut fails = 0;
        if self.initial_solver != InitialSolver::None {
            println!("{}", best_mcm.mcm);
        }
        loop {
            if fails >= max_fails {
                break;
            }
            let original = best_mcm.clone();
            progress.set_description(format!(
                "Choose {amount} - {fails}/{max_fails} - {:.0}",
                best_mcm.log_e
            ));
            if self
                .find_best_check_amount(log_e_cache, progress, amount, &original)
                .and_then(|c| self.update_global_best(best_mcm, &c).then_some(()))
                .is_none()
            {
                fails += 1;
            };
        }
    }

    fn solve_constructive(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        best_mcm: &mut LogeMCM,
    ) -> Bar {
        let mut current_solution =
            MinimallyComplexModel::empty(self.dataset.variables().try_into().unwrap());

        let progress = match self.constructive_strategy {
            ConstructiveStrategy::FrontToBack => self.solve_in_order(
                log_e_cache,
                &mut current_solution,
                (0..self.dataset.variables()).collect::<Vec<usize>>(),
            ),
            ConstructiveStrategy::BackToFront => self.solve_in_order(
                log_e_cache,
                &mut current_solution,
                (0..self.dataset.variables()).rev().collect::<Vec<usize>>(),
            ),
            ConstructiveStrategy::RandomOrder => {
                let mut vector = (0..self.dataset.variables()).collect::<Vec<usize>>();
                vector.shuffle(&mut rand::rng());
                self.solve_in_order(log_e_cache, &mut current_solution, vector)
            }
            ConstructiveStrategy::GreedyOrder => {
                self.solve_greedy_order(log_e_cache, &mut current_solution)
            }
            ConstructiveStrategy::VarianceFirst => {
                self.solve_variance_first(log_e_cache, &mut current_solution)
            }
            ConstructiveStrategy::VarianceLast => {
                self.solve_variance_last(log_e_cache, &mut current_solution)
            }
            _ => panic!("Dont use that one!"),
        };
        *best_mcm = LogeMCM::calculate(current_solution, &self.dataset, log_e_cache);
        progress
    }

    fn solve_in_order(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
        iterator: Vec<usize>,
    ) -> Bar {
        let mut progress = tqdm!(total = iterator.len());
        for var in iterator {
            let vec_rep = current_solution.to_vector();
            let max_icc = vec_rep.iter().max().unwrap() + 1;

            let candidates: Vec<(MinimallyComplexModel, f64)> = (1..=max_icc)
                .map(|c| {
                    let mut new_candidate = vec_rep.clone();
                    new_candidate[var] = c;
                    let mcm = MinimallyComplexModel::from_vector(new_candidate);
                    let log_e = mcm.log_e(&self.dataset, log_e_cache);
                    (mcm, log_e)
                })
                .collect();

            *current_solution = candidates
                .into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap()
                .0;

            progress.set_description(format!(
                "Log E: {:.0}",
                current_solution.log_e(&self.dataset, log_e_cache)
            ));

            progress.update(1).unwrap();
        }

        progress
    }

    fn solve_greedy_order(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
    ) -> Bar {
        let mut progress = tqdm!(total = self.dataset.variables());
        for _ in 0..self.dataset.variables() {
            let vec_rep = current_solution.to_vector();
            let max_icc = vec_rep.iter().max().unwrap() + 1;

            let candidates: Vec<(MinimallyComplexModel, f64)> = vec_rep
                .clone()
                .into_iter()
                .enumerate()
                .filter(|x| x.1 == 0)
                .flat_map(|(x, _)| (1..=max_icc).map(move |c| (x, c)))
                .map(|(x, c)| {
                    let mut new_candidate = vec_rep.clone();
                    new_candidate[x] = c;
                    let mcm = MinimallyComplexModel::from_vector(new_candidate);
                    let log_e = mcm.log_e(&self.dataset, log_e_cache);
                    (mcm, log_e)
                })
                .collect();

            *current_solution = candidates
                .into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap()
                .0;

            progress.update(1).unwrap();
        }

        progress
    }

    fn get_variable_variance(&self) -> Vec<(usize, Vec<usize>, f64)> {
        let mut pairs = vec![];
        for i in 0..self.dataset.variables() {
            let prevalence = self.dataset.state_prevalence(i);

            let count = prevalence.len() as f64;
            let mean = prevalence.iter().sum::<usize>() as f64 / count;
            let variance = prevalence
                .iter()
                .map(|p| (mean - *p as f64).powi(2) / count)
                .sum::<f64>();
            pairs.push((i, prevalence, variance));
        }
        pairs.sort_by(|a, b| a.2.total_cmp(&b.2));
        pairs
    }

    fn solve_variance_first(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
    ) -> Bar {
        let pairs = self.get_variable_variance();
        // pairs.iter().for_each(|p| println!("{p:?}"));

        self.solve_in_order(
            log_e_cache,
            current_solution,
            pairs.iter().map(|x| x.0).collect(),
        )
    }

    fn solve_variance_last(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
    ) -> Bar {
        let pairs = self.get_variable_variance();
        // pairs.iter().for_each(|p| println!("{p:?}"));

        self.solve_in_order(
            log_e_cache,
            current_solution,
            pairs.iter().map(|x| x.0).rev().collect(),
        )
    }
}
fn count_calculations(length: &mut usize, parts_left: usize) {
    for basis in 1..=parts_left {
        for _ in 0..basis {
            *length += 1;
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

impl Solver for GreedySolver {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(GreedySolver {
            dataset: VecDataset::read_from_file(filepath)?,
            continue_after_minimum: false,
            initial_solver: InitialSolver::default(),
            constructive_strategy: ConstructiveStrategy::FrontToBack,
            refinement_sequence: Vec::default(),
        })
    }

    fn solve(&self) -> SolverReport {
        let mut log_e_cache = get_log_e_cache();

        let mut best_mcm = LogeMCM::calculate(
            MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap()),
            &self.dataset,
            &mut log_e_cache,
        );

        let mut progress = match self.initial_solver {
            InitialSolver::Merge => self.solve_merge(&mut log_e_cache, &mut best_mcm, 1),
            InitialSolver::MergeBeam { amount } => {
                self.solve_merge(&mut log_e_cache, &mut best_mcm, amount)
            }
            InitialSolver::Construct => self.solve_constructive(&mut log_e_cache, &mut best_mcm),
            InitialSolver::None => tqdm!(total = 1),
        };

        for refinement in &self.refinement_sequence {
            match refinement {
                Refinements::Local => {
                    self.solve_local(&mut log_e_cache, &mut best_mcm, &mut progress)
                }
                Refinements::ChooseN { n, max_fails } => self.solve_check_amount(
                    *n,
                    *max_fails,
                    &mut log_e_cache,
                    &mut best_mcm,
                    &mut progress,
                ),
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

struct ProductIterator {
    range: usize,
    first: bool,
    current: Vec<usize>,
}

impl ProductIterator {
    fn new(amount: usize, range: usize) -> Self {
        Self {
            range,
            first: true,
            current: vec![0; amount],
        }
    }

    fn len(&self) -> usize {
        (self.range + 1).pow(self.current.len() as u32)
    }
}

impl Iterator for ProductIterator {
    type Item = Vec<usize>;
    fn next(&mut self) -> Option<Self::Item> {
        if !self.first {
            for i in 0..=self.current.len() {
                let x = self.current.get_mut(i)?;
                if *x < self.range {
                    *x += 1;
                    break;
                }
                *x = 0;
            }
        }
        self.first = false;
        Some(self.current.clone())
    }
}
