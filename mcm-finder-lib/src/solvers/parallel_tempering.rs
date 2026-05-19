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

/// How to calculate the temperatures of all of the tempering pools between the
/// minimum (T0) and maximum (Tm) temperature.
pub enum ParallelTemperatureCurve {
    /// Uses linear interpolation between T0 and Tm.
    Linear,
    /// Temperatures calculated using the inverse of the T0 and Tm.
    ///
    /// This one is pretty bad, don't use it.
    InverseLinear,

    Geometric,
    Exponential,
    /// Show you're better than the machines by supplying your own temperatures.
    ///
    /// Requires `self.set_custom_temperatures()` to be set before starting the solve.
    Custom,
}

pub struct ParallelTemperingSearcher {
    dataset: VecDataset,
    starter: AnnealingStarter,
    custom_temperatures: Vec<f64>,
    max_temperature: f64,
    min_temperature: f64,
    temperature_curve: ParallelTemperatureCurve,
    pool_amount: usize,
    steps_per_shuffle: usize,
    shuffles: usize,
    acception_fraction: f64,
}

impl ParallelTemperingSearcher {
    /// Set the starting point for the algorithm
    pub fn set_starter(mut self, starter: AnnealingStarter) -> Self {
        self.starter = starter;
        self
    }

    /// Set the desired temperatures for the tempering pools.
    ///
    /// The amount of temperatures must be at least 2. This option is ignored, unless
    /// `self.set_temperature_curve(ParallelTemperatureCurve::Custom)` is set.
    pub fn set_custom_temperatures(mut self, temperatures: Vec<f64>) -> Self {
        self.custom_temperatures = temperatures;
        self
    }

    /// Set the highest pool temperature.
    ///
    /// If set to `0.0` (default), the solver will calculate it's own max temperature.
    /// Should only be set if you (think you) know better.
    ///
    /// Ignored if `self.set_temperature_curve(ParallelTemperatureCurve::Custom)` is set.
    pub fn set_max_temperature(mut self, temp: f64) -> Self {
        self.max_temperature = temp;
        self
    }

    /// Sets the temperature of the lowers tempering pool. Must be lower than
    /// `self.set_max_temperature()`.
    ///
    /// Ignored if `self.temperature_curve == ParallelTemperatureCurve::Custom` is set.
    pub fn set_min_temperature(mut self, temp: f64) -> Self {
        self.min_temperature = temp;
        self
    }

    /// Sets the way the temperatures of the tempering pools should be distributed.
    ///
    /// Note that for `ParallelTemperingCurve::Custom`, you also need to set
    /// `self.set_custom_temperatures()`
    pub fn set_temperature_curve(mut self, curve: ParallelTemperatureCurve) -> Self {
        self.temperature_curve = curve;
        self
    }

    /// Sets the amount of tempering pools used.
    ///
    /// Default is set to the number of CPU threads, but you could also set this
    /// higher, as the threads wait a lot for access to the log_e_cache.
    ///
    /// Ignored if `self.temperature_curve == ParallelTemperatureCurve::Custom`
    /// is set.
    pub fn set_pools(mut self, pools: usize) -> Self {
        assert!(pools > 1, "Amount of pools should be at least 2.");
        self.pool_amount = pools;
        self
    }

    /// Sets how many steps the tempering pools should temper, before exchanging
    /// MCMs.
    pub fn set_steps_per_shuffle(mut self, steps: usize) -> Self {
        self.steps_per_shuffle = steps;
        self
    }

    /// How many times the MCMs should be exchanged between tempering pools.
    pub fn set_shuffles(mut self, shuffles: usize) -> Self {
        self.shuffles = shuffles;
        self
    }

    /// Sets the target acceptence probability of the dynamically calculated
    /// maximum tempering temperature.
    ///
    /// Ignored when `self.max_temperature()` is set to a non-zero value.
    pub fn set_acception_fraction(mut self, fraction: f64) -> Self {
        self.acception_fraction = fraction;
        self
    }

    fn calculate_inital_temperature(
        &self,
        current: &mut MinimallyComplexModel,
        rng: &mut ThreadRng,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
    ) -> f64 {
        let mut deltas_log_e_regressions: Vec<f64> = vec![];
        let mut progress = par_tqdm!(total = 100, leave = true);
        progress.set_description("Calculating initial temperature");

        while deltas_log_e_regressions.len() < 100 {
            let old_log_e = current.par_log_e(&self.dataset, log_e_cache);
            *current = current.mutate(rng);
            let new_log_e = current.par_log_e(&self.dataset, log_e_cache);

            if new_log_e.total_cmp(&old_log_e).is_lt() {
                deltas_log_e_regressions.push(new_log_e - old_log_e);
                let _ = progress.update(1);
            }
        }

        self.calculate_start_temperature(deltas_log_e_regressions)
    }

    fn calculate_start_temperature(&self, deltas_log_e_regressions: Vec<f64>) -> f64 {
        let average_regression =
            deltas_log_e_regressions.iter().sum::<f64>() / deltas_log_e_regressions.len() as f64;
        average_regression / self.acception_fraction.ln()
    }

    fn linear_temperatures(&self, max_temp: f64) -> Vec<f64> {
        (0..=self.pool_amount)
            .map(|i| {
                self.min_temperature
                    + (max_temp - self.min_temperature) * (i as f64) / (self.pool_amount as f64)
            })
            .collect()
    }

    fn inverse_linear_temperatures(&self, max_temp: f64) -> Vec<f64> {
        (0..=self.pool_amount)
            .map(|i| {
                1.0 / max_temp
                    + (1.0 / self.min_temperature - 1.0 / max_temp) * (i as f64)
                        / (self.pool_amount as f64)
            })
            .map(|t| 1.0 / t)
            .collect()
    }

    fn geometric_temperatures(&self, max_temp: f64) -> Vec<f64> {
        let r = (max_temp / self.min_temperature).powf(1.0 / (self.pool_amount as f64));

        (0..=(self.pool_amount as i32))
            .map(|i| self.min_temperature * r.powi(i))
            .collect()
    }

    fn exponential_temperatures(&self, max_temp: f64) -> Vec<f64> {
        (0..=(self.pool_amount as i32))
            .map(|i| {
                (self.min_temperature.ln()
                    + (max_temp.ln() - self.min_temperature.ln()) * (i as f64)
                        / (self.pool_amount as f64))
                    .exp()
            })
            .collect()
    }

    /// Here we want to use the user provided temperatures when
    /// `ParallelTemperatureCurve::Custom` is provided, and calculate the temperatures
    /// dynamically otherwise. Supplying `self.max_temperature` is optional.
    fn calculate_temperatures(
        &self,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
        starter: &mut MinimallyComplexModel,
    ) -> Vec<f64> {
        let mut temperatures = if let ParallelTemperatureCurve::Custom = self.temperature_curve {
            assert!(
                self.custom_temperatures.len() > 2,
                "At least 2 temperatures are required for custom temperature curve."
            );
            self.custom_temperatures.clone()
        } else {
            let max_temp = if self.max_temperature > 0.0 {
                self.max_temperature
            } else {
                self.calculate_inital_temperature(starter, &mut rand::rng(), log_e_cache)
            };
            match self.temperature_curve {
                ParallelTemperatureCurve::Linear => self.linear_temperatures(max_temp),
                ParallelTemperatureCurve::InverseLinear => {
                    self.inverse_linear_temperatures(max_temp)
                }
                ParallelTemperatureCurve::Geometric => self.geometric_temperatures(max_temp),
                ParallelTemperatureCurve::Exponential => self.exponential_temperatures(max_temp),
                ParallelTemperatureCurve::Custom => {
                    unreachable!("Filtered out in above if statement.")
                }
            }
        };
        temperatures.sort_by(|a, b| b.total_cmp(a));
        temperatures
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
            custom_temperatures: vec![],
            max_temperature: 0.0,
            min_temperature: 0.001,
            temperature_curve: ParallelTemperatureCurve::Linear,
            pool_amount: std::thread::available_parallelism().unwrap().into(),
            steps_per_shuffle: 100,
            shuffles: 200,
            acception_fraction: 0.23,
        })
    }

    fn solve(&self) -> SolverReport {
        let log_e_cache = get_par_log_e_cache();

        let mut starter = match self.starter {
            AnnealingStarter::Single => {
                MinimallyComplexModel::full(NonZero::new(self.dataset.variables()).unwrap())
            }
            AnnealingStarter::Trivial => {
                MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap())
            }
        };

        let temperatures = self.calculate_temperatures(&log_e_cache, &mut starter);
        println!("{:?}", temperatures);

        let mut pools: Vec<_> = temperatures
            .iter()
            .map(|t| TemperPool::new(starter.clone(), *t, &self.dataset))
            .collect();

        let bar = Arc::new(Mutex::new(par_tqdm!(
            total = self.shuffles * self.steps_per_shuffle * temperatures.len()
        )));

        for i in 0..self.shuffles {
            pools
                .par_iter_mut()
                .for_each(|p| p.step_n(&self.dataset, self.steps_per_shuffle, &log_e_cache, &bar));
            let best_pool = pools
                .iter()
                .max_by(|a, b| a.best_log_e.total_cmp(&b.best_log_e))
                .unwrap();
            bar.lock()
                .unwrap()
                .set_description(format!("Best Log(E): {:.0}", best_pool.best_log_e));

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
            .max_by(|a, b| a.best_log_e.total_cmp(&b.best_log_e))
            .unwrap();

        SolverReport::new(
            best_pool.best_mcm.clone(),
            best_pool.best_log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{:.0}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}

pub struct TemperPool {
    mcm: MinimallyComplexModel,
    log_e: f64,
    temperature: f64,
    best_mcm: MinimallyComplexModel,
    best_log_e: f64,
}

impl TemperPool {
    pub fn new(mcm: MinimallyComplexModel, temperature: f64, dataset: &VecDataset) -> Self {
        let log_e = mcm.log_e(dataset, &mut None);
        Self {
            best_mcm: mcm.clone(),
            best_log_e: log_e,
            mcm,
            log_e,
            temperature,
        }
    }

    /// Perform a Metropolis step on the held `MinimallyComplexModel`.
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
            if new_log_e.total_cmp(&self.best_log_e).is_gt() {
                self.best_log_e = new_log_e;
                self.best_mcm = new_mcm.clone();
            }

            self.mcm = new_mcm;
            self.log_e = new_log_e;
        }
    }

    /// Perform n Metropolis steps. Also updates the progress bar.
    fn step_n(
        &mut self,
        dataset: &VecDataset,
        n: usize,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
        bar: &Arc<Mutex<Bar>>,
    ) {
        for _ in 0..n {
            rayon::join(
                || {
                    let mut rng = rand::rng();
                    self.step(dataset, &mut rng, log_e_cache);
                },
                || bar.lock().unwrap().update(1).unwrap(),
            );
        }
    }

    /// Stochatically swaps the `MinimallyComplesModel`s between this and the supplied
    /// `TemperPool`.
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
