use std::{
    collections::HashMap,
    mem,
    num::NonZero,
    path::Path,
    sync::{Arc, Mutex},
};

use dashmap::DashMap;
use fixedbitset::FixedBitSet;
use kdam::{Bar, BarExt, par_tqdm};
use rand::{RngExt, rngs::ThreadRng};
use rayon::prelude::*;

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::{
        AnnealingStarter, get_par_log_e_cache,
        solvers_base::{Solver, SolverReport},
    },
};

pub struct ParallelTemperingSearcher {
    dataset: VecDataset,
    starter: AnnealingStarter,
    temperatures: Vec<f64>,
    steps_per_shuffle: usize,
    shuffles: usize,
}

impl ParallelTemperingSearcher {
    pub fn set_starter(mut self, starter: AnnealingStarter) -> Self {
        self.starter = starter;
        self
    }

    pub fn set_temperatures(mut self, temperatures: Vec<f64>) -> Self {
        self.temperatures = temperatures;
        self
    }
}

impl Solver for ParallelTemperingSearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(ParallelTemperingSearcher {
            dataset: VecDataset::read_from_file(filepath)?,
            starter: AnnealingStarter::default(),
            temperatures: vec![],
            steps_per_shuffle: 100,
            shuffles: 100,
        })
    }

    fn solve(&self) -> SolverReport {
        let starter = match self.starter {
            AnnealingStarter::Single => {
                MinimallyComplexModel::full(NonZero::new(self.dataset.variables()).unwrap())
            }
            AnnealingStarter::Trivial => {
                MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap())
            }
        };

        let mut pools: Vec<_> = self
            .temperatures
            .iter()
            .map(|t| TemperPool::new(starter.clone(), *t, &self.dataset))
            .collect();

        let log_e_cache = get_par_log_e_cache();
        let bar = Arc::new(Mutex::new(par_tqdm!(
            total = self.shuffles * self.steps_per_shuffle * self.temperatures.len()
        )));

        for i in 0..self.shuffles {
            pools.par_iter_mut().for_each_init(rand::rng, |rng, p| {
                p.step_n(
                    &self.dataset,
                    rng,
                    self.steps_per_shuffle,
                    &log_e_cache,
                    &bar,
                )
            });
            if i % 2 == 0 {
                pools
                    .par_chunks_exact_mut(2)
                    .for_each_init(rand::rng, |rng, p| {
                        let (p_one, p_two) = p.split_at_mut(1);
                        p_one[0].swap(&mut p_two[0], rng);
                    });
            }
        }

        let best_pool = pools
            .iter()
            .max_by(|a, b| a.log_e.total_cmp(&b.log_e))
            .unwrap();

        SolverReport::new(
            best_pool.mcm.clone(),
            best_pool.log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}

pub struct TemperPool {
    mcm: MinimallyComplexModel,
    log_e: f64,
    temperature: f64,
}

impl TemperPool {
    pub fn new(mcm: MinimallyComplexModel, temperature: f64, dataset: &VecDataset) -> Self {
        let log_e = mcm.log_e(dataset, &mut None);
        Self {
            mcm,
            log_e,
            temperature,
        }
    }

    fn step(
        &mut self,
        dataset: &VecDataset,
        rng: &mut ThreadRng,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
    ) {
        let new_mcm = self.mcm.mutate(rng);
        let new_log_e = new_mcm.par_log_e(dataset, log_e_cache);
        if new_log_e.total_cmp(&self.log_e).is_gt()
            || rng.random_bool(((new_log_e - self.log_e) / self.temperature).exp())
        {
            self.mcm = new_mcm;
            self.log_e = new_log_e;
        }
    }

    fn step_n(
        &mut self,
        dataset: &VecDataset,
        rng: &mut ThreadRng,
        n: usize,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
        bar: &Arc<Mutex<Bar>>,
    ) {
        for _ in 0..n {
            self.step(dataset, rng, log_e_cache);
            let _ = bar.lock().unwrap().update(1);
        }
    }

    fn swap(&mut self, other: &mut TemperPool, rng: &mut ThreadRng) {
        let probability = ((1.0 / other.temperature - 1.0 / self.temperature)
            * (other.log_e - self.log_e))
            .exp()
            .min(1.0);
        if rng.random_bool(probability) {
            mem::swap(&mut self.mcm, &mut other.mcm);
            mem::swap(&mut self.log_e, &mut other.log_e);
        }
    }
}
