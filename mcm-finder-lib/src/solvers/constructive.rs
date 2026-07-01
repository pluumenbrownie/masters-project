use std::{
    collections::{HashMap, VecDeque},
    marker::Sized,
    ops::{Add, Div},
    path::Path,
};

use fixedbitset::FixedBitSet;
use kdam::{BarExt, tqdm};
use rand::{RngExt, rngs::ThreadRng, seq::SliceRandom};

use crate::{
    dataset::{Dataset, simple::VecDataset},
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
    /// Iteratively constructs MCM in a random order, while building a tabu list
    /// for good and bad variable combinations.
    Tabu {
        steps: usize,
    },
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
            ConstructiveStrategy::Tabu { steps } => {
                self.solve_tabu(&mut log_e_cache, &mut current_solution, steps)
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

    fn solve_tabu(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
        steps: usize,
    ) {
        let rng = &mut rand::rng();
        let mut order_vector = (0..self.dataset.variables()).collect::<Vec<usize>>();

        let mut tabu_list = vec![vec![0i8; self.dataset.variables()]; self.dataset.variables()];
        let mut past_log_e = LogEQueue::default();

        let mut bar = tqdm!(total = steps, position = 1);
        for _step in 0..steps {
            *current_solution =
                MinimallyComplexModel::empty(self.dataset.variables().try_into().unwrap());
            order_vector.shuffle(rng);
            self.solve_in_order_tabu(
                log_e_cache,
                current_solution,
                &order_vector,
                &tabu_list,
                rng,
            );
            let improvement =
                past_log_e.improvement(current_solution.log_e(&self.dataset, log_e_cache));
            if let Some(imp) = improvement {
                update_tabu(&mut tabu_list, current_solution, imp);
            }
            // println!("{:?}", tabu_list);
            bar.set_description(format!(
                "Log E: {:.0}",
                current_solution.log_e(&self.dataset, log_e_cache)
            ));
            bar.update(1).unwrap();
        }
    }

    fn solve_in_order_tabu(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        current_solution: &mut MinimallyComplexModel,
        iterator: &[usize],
        tabu_list: &[Vec<i8>],
        rng: &mut ThreadRng,
    ) {
        for var in tqdm!(iterator.iter()) {
            let vec_rep = current_solution.to_vector();
            let max_icc = vec_rep.iter().max().unwrap() + 1;
            let var_amount = vec_rep.iter().filter(|x| **x != 0).count() + 1;

            let mut candidates: Vec<(MinimallyComplexModel, f64, f64)> = (1..=max_icc)
                .map(|new_icc| {
                    let mut new_candidate = vec_rep.clone();
                    new_candidate[*var] = new_icc;
                    let mcm = MinimallyComplexModel::from_vector(new_candidate);
                    let log_e = mcm.log_e(&self.dataset, log_e_cache);
                    let tabu_count: f64 = vec_rep
                        .iter()
                        .enumerate()
                        .filter_map(|(v, icc)| {
                            if icc == &new_icc {
                                Some(tabu_list[*var][v] as f64)
                            } else {
                                None
                            }
                        })
                        .sum();
                    let goodness = 1.0 + (tabu_count / 2.0).tanh().min(0.0);
                    (mcm, log_e, goodness)
                })
                .collect();

            candidates.sort_by(|(_, a, _), (_, b, _)| b.total_cmp(a));
            let goodness_vec: Vec<f64> = candidates.iter().map(|x| x.2).collect();
            if var_amount == 61 {
                println!("{:?}", goodness_vec);
            }

            let mut counter = 0;
            *current_solution = loop {
                if counter >= max_icc {
                    counter = 0
                };
                let goodness = candidates[counter].2;
                if rng.random_bool(goodness) {
                    break candidates[counter].0.to_owned();
                }
                counter += 1;
            };
        }
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

fn update_tabu(tabu_list: &mut [Vec<i8>], current_solution: &MinimallyComplexModel, imp: bool) {
    let matrix = current_solution.to_matrix();
    for (tabu_row, mcm_row) in tabu_list.iter_mut().zip(matrix.iter()) {
        for c in mcm_row.ones() {
            tabu_row[c] = tabu_row[c].add(if imp { 1 } else { -1 }).clamp(-3, 3);
        }
    }
}

#[derive(Debug, Default)]
struct LogEQueue {
    data: [Option<f64>; 16],
    counter: u8,
}

impl LogEQueue {
    fn add(&mut self, value: f64) {
        self.data[(self.counter % 16) as usize] = Some(value);
        self.counter = self.counter.wrapping_add(1);
    }

    fn avg(&self) -> Option<f64> {
        self.data[0]?;
        let total: f64 = self.data.iter().filter_map(|x| *x).sum();
        Some(total / self.data.iter().filter(|x| x.is_some()).count() as f64)
    }

    /// Adds the new log E, returns if the new value is better or worse than average.
    ///
    /// Returns `None` when the new value is the first one added.
    #[must_use]
    fn improvement(&mut self, value: f64) -> Option<bool> {
        let avg = self.avg();
        self.add(value);
        avg.map(|a| value.total_cmp(&a).is_gt())
    }
}
