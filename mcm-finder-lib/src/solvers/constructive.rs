use std::{collections::HashMap, marker::Sized, path::Path};

use fixedbitset::FixedBitSet;
use kdam::tqdm;
use rand::seq::SliceRandom;

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
pub enum ConstructiveStrategy {
    /// Start with variable 0, then place variables incrementally.
    FrontToBack,
    /// Start with last variable, then place variables decrementally.
    BackToFront,
    /// Place variables in random order.
    RandomOrder,
    /// Greedily find the best variable to place in each step.
    GreedyOrder,
    VarianceFirst,
    VarianceLast,
}

#[derive(Clone)]
pub struct ConstructiveSolver {
    dataset: VecDataset,
    strategy: ConstructiveStrategy,
}

impl Solver for ConstructiveSolver {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(ConstructiveSolver {
            dataset: VecDataset::read_from_file(filepath)?,
            strategy: ConstructiveStrategy::FrontToBack,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut log_e_cache = get_log_e_cache();

        let mut current_solution =
            MinimallyComplexModel::empty(self.dataset.variables().try_into().unwrap());

        match self.strategy {
            ConstructiveStrategy::FrontToBack => self.solve_in_order(
                &mut log_e_cache,
                &mut current_solution,
                (0..self.dataset.variables()).collect::<Vec<usize>>(),
            ),
            ConstructiveStrategy::BackToFront => self.solve_in_order(
                &mut log_e_cache,
                &mut current_solution,
                (0..self.dataset.variables()).rev().collect::<Vec<usize>>(),
            ),
            ConstructiveStrategy::RandomOrder => {
                let mut vector = (0..self.dataset.variables()).collect::<Vec<usize>>();
                vector.shuffle(&mut rand::rng());
                self.solve_in_order(&mut log_e_cache, &mut current_solution, vector)
            }
            ConstructiveStrategy::GreedyOrder => {
                self.solve_greedy_order(&mut log_e_cache, &mut current_solution)
            }
            ConstructiveStrategy::VarianceFirst => {
                self.solve_variance_first(&mut log_e_cache, &mut current_solution)
            }
            ConstructiveStrategy::VarianceLast => {
                self.solve_variance_last(&mut log_e_cache, &mut current_solution)
            }
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

impl ConstructiveSolver {
    pub fn set_strategy(mut self, strategy: ConstructiveStrategy) -> ConstructiveSolver {
        self.strategy = strategy;
        self
    }

    fn solve_in_order(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
        iterator: Vec<usize>,
    ) {
        for var in tqdm!(iterator.iter()) {
            let vec_rep = current_solution.to_vector();
            let max_icc = vec_rep.iter().max().unwrap() + 1;

            let candidates: Vec<(MinimallyComplexModel, f64)> = (1..=max_icc)
                .map(|c| {
                    let mut new_candidate = vec_rep.clone();
                    new_candidate[*var] = c;
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
        }
    }

    fn solve_greedy_order(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
    ) {
        for _ in tqdm!(0..self.dataset.variables()) {
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
        }
    }

    fn solve_variance_first(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
    ) {
        let pairs = self.get_variable_variance();
        // pairs.iter().for_each(|p| println!("{p:?}"));

        self.solve_in_order(
            log_e_cache,
            current_solution,
            pairs.iter().map(|x| x.0).collect(),
        );
    }

    fn solve_variance_last(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
    ) {
        let pairs = self.get_variable_variance();
        // pairs.iter().for_each(|p| println!("{p:?}"));

        self.solve_in_order(
            log_e_cache,
            current_solution,
            pairs.iter().map(|x| x.0).rev().collect(),
        );
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
}
