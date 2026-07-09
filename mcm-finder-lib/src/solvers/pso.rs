use std::{cell::RefCell, collections::HashMap, marker::Sized, path::Path, rc::Rc};

use fixedbitset::FixedBitSet;
use kdam::{BarExt, tqdm};
use rand::{RngExt, rngs::ThreadRng};

use crate::{
    dataset::{Dataset, simple::VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::{
        get_log_e_cache,
        solvers_base::{Solver, SolverReport},
    },
};

#[derive(Clone, Copy)]
pub enum PsoStrategy {
    Rounding,
    ConstructiveOrder,
}

#[derive(Clone)]
pub struct PsoSolver {
    dataset: VecDataset,
    swarm_size: usize,
    starter_iccs: Option<usize>,
    strategy: PsoStrategy,
    momentum_weight: f64,
    cognitive_weight: f64,
    social_weight: f64,
    steps: usize,
    particle_max_weight: f64,
}

impl Solver for PsoSolver {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(PsoSolver {
            dataset: VecDataset::read_from_file(filepath)?,
            swarm_size: 64,
            starter_iccs: None,
            strategy: PsoStrategy::Rounding,
            momentum_weight: 0.8,
            cognitive_weight: 1.4,
            social_weight: 1.45,
            steps: 1000,
            particle_max_weight: 20.0,
        })
    }

    fn solve(&self) -> SolverReport {
        let log_e_cache = Rc::new(RefCell::new(get_log_e_cache()));
        let global_best = Rc::new(RefCell::new(None));

        match self.strategy {
            PsoStrategy::Rounding => self.solve_rounding(&log_e_cache, &global_best),
            PsoStrategy::ConstructiveOrder => {
                self.solve_constuctive_order(&log_e_cache, &global_best)
            }
        }

        let best_mcm = float_to_mcm(&global_best.borrow().clone().unwrap());
        let log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache.borrow_mut());
        SolverReport::new(
            best_mcm,
            log_e,
            HashMap::from([
                (
                    "Unique ICCs covered".into(),
                    format!("{}", log_e_cache.borrow().as_ref().unwrap().len()),
                ),
                (
                    "Best Particle X".into(),
                    format!("{:?}", global_best.borrow().as_ref().unwrap()),
                ),
            ]),
        )
    }
}

impl PsoSolver {
    pub fn set_momentum_weight(mut self, weight: f64) -> Self {
        self.momentum_weight = weight;
        self
    }

    pub fn set_cognitive_weight(mut self, weight: f64) -> Self {
        self.cognitive_weight = weight;
        self
    }

    pub fn set_social_weight(mut self, weight: f64) -> Self {
        self.social_weight = weight;
        self
    }
    pub fn set_strategy(mut self, strategy: PsoStrategy) -> Self {
        self.strategy = strategy;
        self
    }

    pub fn set_swarm_size(mut self, swarm_size: usize) -> Self {
        self.swarm_size = swarm_size;
        self
    }

    pub fn set_steps(mut self, steps: usize) -> Self {
        self.steps = steps;
        self
    }

    fn get_mcm_range(&self) -> usize {
        self.starter_iccs
            .unwrap_or_else(|| self.dataset.variables())
    }

    fn solve_rounding(
        &self,
        log_e_cache: &Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
        global_best: &Rc<RefCell<Option<Vec<f64>>>>,
    ) {
        let mut rng = rand::rng();
        let mut swarm = vec![];
        for _ in 0..self.swarm_size {
            swarm.push(Particle::new(
                self,
                global_best.clone(),
                log_e_cache.clone(),
                &mut rng,
            ));
        }

        let mut progress = tqdm!(total = self.steps);
        for _ in 0..self.steps {
            swarm.iter_mut().for_each(|p| p.evaluate(&self.dataset));
            swarm.iter_mut().for_each(|p| p.update(&mut rng));

            progress.set_description(format!(
                "Log(E) = {:.0}",
                float_to_mcm(&global_best.borrow().clone().unwrap())
                    .log_e(&self.dataset, &mut log_e_cache.borrow_mut())
            ));
            let _ = progress.update(1);
        }
    }

    fn solve_constuctive_order(
        &self,
        log_e_cache: &Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
        global_best: &Rc<RefCell<Option<Vec<f64>>>>,
    ) {
        let mut rng = rand::rng();
        let mut swarm = vec![];
        for _ in 0..self.swarm_size {
            swarm.push(Particle::new(
                self,
                global_best.clone(),
                log_e_cache.clone(),
                &mut rng,
            ));
        }

        let mut progress = tqdm!(total = self.steps);
        for _ in 0..self.steps {
            swarm.iter_mut().for_each(|p| p.evaluate(&self.dataset));
            swarm.iter_mut().for_each(|p| p.update(&mut rng));

            progress.set_description(format!(
                "Log(E) = {:.0}",
                float_to_mcm(&global_best.borrow().clone().unwrap())
                    .log_e(&self.dataset, &mut log_e_cache.borrow_mut())
            ));
            let _ = progress.update(1);
        }
    }
}

#[derive(Clone)]
struct Weights {
    momentum: f64,
    cognitive: f64,
    social: f64,
}

#[derive(Clone)]
struct Particle {
    weights: Weights,
    x: Vec<f64>,
    v: Option<Vec<f64>>,
    best_x: Option<Vec<f64>>,
    global_best_x: Rc<RefCell<Option<Vec<f64>>>>,
    cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
    strategy: PsoStrategy,
}

impl Particle {
    fn new(
        solver: &PsoSolver,
        global_best: Rc<RefCell<Option<Vec<f64>>>>,
        cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
        rng: &mut ThreadRng,
    ) -> Particle {
        let numbers = Vec::from_iter((0..solver.dataset.variables()).map(|_| {
            rng.random_range(solver.particle_max_weight..(solver.particle_max_weight + 5.0))
        }));
        // let numbers = vec![5.0f64; solver.dataset.variables()];
        let mut particle = Particle {
            weights: Weights {
                momentum: solver.momentum_weight,
                cognitive: solver.cognitive_weight,
                social: solver.social_weight,
            },
            x: numbers,
            v: None,
            best_x: None,
            global_best_x: global_best,
            cache,
            strategy: solver.strategy,
        };
        particle.evaluate(&solver.dataset);
        particle
    }

    fn evaluate(&mut self, dataset: &VecDataset) {
        let mcm = self.generate_mcm(dataset);

        let log_e = mcm.log_e(dataset, &mut self.cache.borrow_mut());

        if self.best_x.as_ref().is_none_or(|b| {
            log_e
                .total_cmp(&float_to_mcm(b).log_e(dataset, &mut self.cache.borrow_mut()))
                .is_gt()
        }) {
            self.best_x = Some(self.x.clone());

            if self.global_best_x.borrow().as_ref().is_none_or(|b| {
                log_e
                    .total_cmp(&float_to_mcm(b).log_e(dataset, &mut self.cache.borrow_mut()))
                    .is_gt()
            }) {
                self.global_best_x.replace(Some(self.x.clone()));
            }
        }
    }

    fn generate_mcm(&mut self, dataset: &VecDataset) -> MinimallyComplexModel {
        match self.strategy {
            PsoStrategy::Rounding => float_to_mcm(&self.x),
            PsoStrategy::ConstructiveOrder => self.constructive_mcm(dataset),
        }
    }

    fn update(&mut self, rng: &mut ThreadRng) {
        let velocity = self.get_velocity(rng);

        self.x.iter_mut().zip(velocity.iter()).for_each(|(x, v)| {
            *x += v;
            *x = x.clamp(0.51, f64::MAX);
        });
        self.v = Some(velocity);
    }

    fn get_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        let momentum = self.get_momentum_velocity(rng);
        let cognitive = self.get_cognitive_velocity(rng);
        let social = self.get_social_velocity(rng);

        let velocity: Vec<f64> = momentum
            .iter()
            .zip(cognitive.iter())
            .zip(social.iter())
            .map(|((a, b), c)| a + b + c)
            .collect();
        velocity
    }

    fn get_momentum_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        let mut v = self.v.clone().unwrap_or_else(|| {
            let mut v = vec![0.0f64; self.x.len()];
            v.iter_mut()
                .for_each(|value| *value = rng.random_range(-5.0f64..5.0));
            v
        });
        v.iter_mut().for_each(|v| *v *= self.weights.momentum);
        v
    }

    fn get_cognitive_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        self.best_x
            .as_deref()
            .unwrap()
            .iter()
            .zip(self.x.iter())
            .map(|(b, x)| (b - x) * rng.random_range(0.0f64..1.0) * self.weights.cognitive)
            .collect()
    }

    fn get_social_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        self.global_best_x
            .borrow()
            .as_deref()
            .unwrap()
            .iter()
            .zip(self.x.iter())
            .map(|(b, x)| (b - x) * rng.random_range(0.0f64..1.0) * self.weights.social)
            .collect()
    }

    fn constructive_mcm(&self, dataset: &VecDataset) -> MinimallyComplexModel {
        let mut current_solution = MinimallyComplexModel::empty(self.x.len().try_into().unwrap());

        let mut iterator: Vec<_> = self.x.iter().enumerate().collect();
        iterator.sort_by(|a, b| a.1.total_cmp(b.1));

        for var in iterator.iter().map(|x| x.0) {
            let vec_rep = current_solution.to_vector();
            let max_icc = vec_rep.iter().max().unwrap() + 1;

            let candidates: Vec<(MinimallyComplexModel, f64)> = (1..=max_icc)
                .map(|c| {
                    let mut new_candidate = vec_rep.clone();
                    new_candidate[var] = c;
                    let mcm = MinimallyComplexModel::from_vector(new_candidate);
                    let log_e = mcm.log_e(dataset, &mut self.cache.borrow_mut());
                    (mcm, log_e)
                })
                .collect();

            current_solution = candidates
                .into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap()
                .0;
        }

        current_solution
    }
}

fn float_to_mcm(x: &[f64]) -> MinimallyComplexModel {
    MinimallyComplexModel::from_vector(x.iter().map(|x| x.round() as usize).collect())
}
