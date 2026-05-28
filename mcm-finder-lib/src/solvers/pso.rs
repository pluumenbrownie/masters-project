use std::{cell::RefCell, collections::HashMap, marker::Sized, path::Path, rc::Rc};

use fixedbitset::FixedBitSet;
use kdam::tqdm;
use rand::{RngExt, rngs::ThreadRng};

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
pub struct PsoSolver {
    dataset: VecDataset,
    swarm_size: usize,
    starter_iccs: Option<usize>,
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
        })
    }

    fn solve(&self) -> SolverReport {
        let log_e_cache = Rc::new(RefCell::new(get_log_e_cache()));
        let mut rng = rand::rng();

        let mut swarm = vec![];
        let global_best = Rc::new(RefCell::new(None));
        for _ in 0..self.swarm_size {
            swarm.push(Particle::new(
                self.dataset.variables(),
                global_best.clone(),
                self.get_mcm_range() as f64,
                log_e_cache.clone(),
                &self.dataset,
                &mut rng,
            ));
        }

        for _ in tqdm!(0usize..10_000) {
            swarm.iter_mut().for_each(|p| p.evaluate(&self.dataset));
            swarm.iter_mut().for_each(|p| p.update(&mut rng));
        }

        let best_mcm = float_to_mcm(&global_best.borrow().clone().unwrap());
        let log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache.borrow_mut());
        SolverReport::new(
            best_mcm,
            log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.borrow().as_ref().unwrap().len()),
            )]),
        )
    }
}

impl PsoSolver {
    fn get_mcm_range(&self) -> usize {
        self.starter_iccs
            .unwrap_or_else(|| self.dataset.variables())
    }
}

#[derive(Clone)]
struct Particle {
    momentum_weight: f64,
    cognitive_weight: f64,
    social_weight: f64,
    x: Vec<f64>,
    v: Option<Vec<f64>>,
    best_x: Option<Vec<f64>>,
    global_best_x: Rc<RefCell<Option<Vec<f64>>>>,
    cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
}

impl Particle {
    fn new(
        variables: usize,
        global_best: Rc<RefCell<Option<Vec<f64>>>>,
        max_value: f64,
        cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
        dataset: &VecDataset,
        rng: &mut ThreadRng,
    ) -> Particle {
        let numbers = Vec::from_iter((0..variables).map(|_| rng.random_range(1.0..max_value)));
        let mut particle = Particle {
            momentum_weight: 0.5,
            cognitive_weight: 0.50,
            social_weight: 0.50,
            x: numbers,
            v: None,
            best_x: None,
            global_best_x: global_best,
            cache,
        };
        particle.evaluate(dataset);
        particle
    }

    fn evaluate(&mut self, dataset: &VecDataset) {
        let mcm = float_to_mcm(&self.x);

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
        // if log_e
        //     .total_cmp(&best_mcm.log_e(dataset, &mut self.cache.borrow_mut()))
        //     .is_gt()
        // {
        //     self.best_x = self.x.clone();

        //     let global_mcm = self.global_best_to_mcm();
        //     match global_mcm {
        //         None => {
        //             self.global_best_x.replace(Some(self.x.clone()));
        //         }
        //         Some(g_mcm) => {
        //             if log_e
        //                 .total_cmp(&g_mcm.log_e(dataset, &mut self.cache.borrow_mut()))
        //                 .is_gt()
        //             {
        //                 self.global_best_x.replace(Some(self.x.clone()));
        //             }
        //         }
        //     }
        // }
    }

    fn update(&mut self, rng: &mut ThreadRng) {
        let velocity = self.get_velocity(rng);

        self.x.iter_mut().zip(velocity.iter()).for_each(|(x, v)| {
            *x += v;
            *x = x.clamp(1.0, f64::MAX);
        });
        self.v = Some(velocity);
    }

    fn get_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        let momentum = self.get_momentum_velocity();
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

    fn get_momentum_velocity(&self) -> Vec<f64> {
        let mut v = self.v.clone().unwrap_or(vec![0.0f64; self.x.len()]);
        v.iter_mut().for_each(|v| *v *= self.momentum_weight);
        v
    }

    fn get_cognitive_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        self.best_x
            .as_deref()
            .unwrap()
            .iter()
            .zip(self.x.iter())
            .map(|(b, x)| (b - x) * rng.random_range(0.0f64..1.0) * self.cognitive_weight)
            .collect()
    }

    fn get_social_velocity(&self, rng: &mut ThreadRng) -> Vec<f64> {
        self.global_best_x
            .borrow()
            .as_deref()
            .unwrap()
            .iter()
            .zip(self.x.iter())
            .map(|(b, x)| (b - x) * rng.random_range(0.0f64..1.0) * self.social_weight)
            .collect()
    }

    fn get_log_e(&self, dataset: &VecDataset) -> f64 {
        let mcm = float_to_mcm(&self.x);
        mcm.log_e(dataset, &mut self.cache.borrow_mut())
    }
}

fn float_to_mcm(x: &[f64]) -> MinimallyComplexModel {
    MinimallyComplexModel::from_vector(x.iter().map(|x| x.round() as usize).collect())
}
