use std::{cell::RefCell, collections::HashMap, num::NonZeroUsize, rc::Rc};

use fixedbitset::FixedBitSet;
use kdam::{BarExt, tqdm};
use rand::{RngExt, rngs::ThreadRng};

use crate::{
    dataset::{Dataset, simple::VecDataset},
    logger::SolverEventSender,
    mcm::MinimallyComplexModel,
    solvers::{Solver, SolverReport, get_log_e_cache},
};

pub struct SpsoSolver {
    dataset: VecDataset,
    swarm_size: usize,
    steps: usize,
    weights: Weights,
    sender: Option<SolverEventSender>,
}

impl Solver for SpsoSolver {
    fn from_file(filepath: &std::path::Path) -> Result<Self, crate::mcm_error::MCMError>
    where
        Self: Sized,
    {
        Ok(SpsoSolver {
            dataset: VecDataset::read_from_file(filepath)?,
            swarm_size: 64,
            steps: 100000,
            weights: Weights { weight: 0.8 },
            sender: None,
        })
    }

    fn solve(&self) -> super::SolverReport {
        let log_e_cache = Rc::new(RefCell::new(get_log_e_cache()));
        let global_best = Rc::new(RefCell::new(None));
        let mut rng = rand::rng();

        let mut swarm = [Particle::init(
            self.dataset.variables().try_into().unwrap(),
            self.weights,
            global_best.clone(),
            log_e_cache.clone(),
            &mut rng,
        )];

        let mut progress = tqdm!(total = self.steps);
        for _ in 0..self.steps {
            swarm.iter_mut().for_each(|p| p.evaluate(&self.dataset));
            swarm.iter_mut().for_each(|p| p.update(&mut rng));

            progress.set_description(format!(
                "Log(E) = {:.0}",
                global_best
                    .borrow()
                    .clone()
                    .unwrap()
                    .log_e(&self.dataset, &mut log_e_cache.borrow_mut())
            ));
            let _ = progress.update(1);
        }

        let best_mcm = &global_best.borrow().clone().unwrap();
        let log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache.borrow_mut());
        SolverReport::new(
            best_mcm.clone(),
            log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.borrow().as_ref().unwrap().len()),
            )]),
        )
    }

    fn get_sender(&self) -> Option<&SolverEventSender> {
        self.sender.as_ref()
    }
}

#[derive(Debug, Clone, Copy)]
struct Weights {
    weight: f64,
}

struct Particle {
    weights: Weights,
    x: MinimallyComplexModel,
    random_x: MinimallyComplexModel,
    best_x: Option<MinimallyComplexModel>,
    global_best_x: Rc<RefCell<Option<MinimallyComplexModel>>>,
    cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
}

impl Particle {
    fn new(
        weights: Weights,
        x: MinimallyComplexModel,
        random_x: MinimallyComplexModel,
        best_x: Option<MinimallyComplexModel>,
        global_best_x: Rc<RefCell<Option<MinimallyComplexModel>>>,
        cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
    ) -> Particle {
        Particle {
            weights,
            x,
            random_x,
            best_x,
            global_best_x,
            cache,
        }
    }

    fn init(
        variables: NonZeroUsize,
        weights: Weights,
        global_best_x: Rc<RefCell<Option<MinimallyComplexModel>>>,
        cache: Rc<RefCell<Option<HashMap<FixedBitSet, f64>>>>,
        rng: &mut ThreadRng,
    ) -> Particle {
        let mut mcm = MinimallyComplexModel::full(variables);
        let mut random_mcm = MinimallyComplexModel::full(variables);
        for _ in 0..15 {
            (mcm, _) = mcm.mutate(rng);
            (random_mcm, _) = random_mcm.mutate(rng);
        }
        Particle::new(weights, mcm, random_mcm, None, global_best_x, cache)
    }

    fn evaluate(&mut self, dataset: &VecDataset) {
        let log_e = self.x.log_e(dataset, &mut self.cache.borrow_mut());

        if self.best_x.as_ref().is_none_or(|b| {
            log_e
                .total_cmp(&b.log_e(dataset, &mut self.cache.borrow_mut()))
                .is_gt()
        }) {
            self.best_x = Some(self.x.clone());

            if self.global_best_x.borrow().as_ref().is_none_or(|b| {
                log_e
                    .total_cmp(&b.log_e(dataset, &mut self.cache.borrow_mut()))
                    .is_gt()
            }) {
                self.global_best_x.replace(Some(self.x.clone()));
            }
        }
    }

    fn update(&mut self, rng: &mut ThreadRng) {
        let distance_random = self.x.jakkard_difference(&self.random_x);
        let distance_best = self.x.jakkard_difference(self.best_x.as_ref().unwrap());
        let distance_global = self
            .x
            .jakkard_difference(&self.global_best_x.borrow().clone().unwrap());
        let mut counter = 0usize;
        self.x = loop {
            let (candidate, _) = self.x.mutate(rng);
            let random_closer =
                distance_random - candidate.jakkard_difference(&self.random_x) < 0.0;
            let best_closer =
                distance_best - candidate.jakkard_difference(self.best_x.as_ref().unwrap()) < 0.0;
            let globar_closer = distance_global
                - candidate.jakkard_difference(&self.global_best_x.borrow().clone().unwrap())
                < 0.0;
            if random_closer && best_closer && globar_closer || counter > 10000 {
                break candidate;
            }
            counter += 1;
        };
        counter = 0;
        if rng.random_bool(self.weights.weight) {
            let distance_best = self.x.jakkard_difference(self.best_x.as_ref().unwrap());
            let distance_global = self
                .x
                .jakkard_difference(&self.global_best_x.borrow().clone().unwrap());
            self.random_x = loop {
                let (candidate, _) = self.random_x.mutate(rng);
                if distance_best - candidate.jakkard_difference(self.best_x.as_ref().unwrap()) < 0.0
                    && distance_global
                        - candidate
                            .jakkard_difference(&self.global_best_x.borrow().clone().unwrap())
                        < 0.0
                    || counter > 10000
                {
                    break candidate;
                }
                counter += 1;
            }
        }
    }
}
