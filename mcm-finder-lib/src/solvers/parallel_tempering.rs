use std::{
    collections::HashMap,
    mem,
    num::NonZero,
    path::Path,
    sync::{Arc, Mutex, mpsc::Sender},
};

use annolog::{Collector, CollectorEvent};
use dashmap::DashMap;
use fixedbitset::FixedBitSet;
use kdam::{Bar, BarExt, par_tqdm};
use rand::{
    RngExt,
    rngs::ThreadRng,
    seq::{IteratorRandom, SliceRandom},
};
use rayon::prelude::*;
use serde::Serialize;

use crate::{
    dataset::{Dataset, simple::VecDataset},
    logger::SolverEvent,
    mcm::{MCMData, MinimallyComplexModel, MutationEvent, ParLogEResult},
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

pub struct ParallelTemperingSolver {
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
    arbitrary_swaps: bool,
    count_non_cached: bool,
    sender: Option<Sender<CollectorEvent<SolverEvent>>>,
}

impl ParallelTemperingSolver {
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
        assert!(temp >= 0.0, "Maximum temperature cannot be negative.");
        assert!(
            temp >= self.min_temperature || temp == 0.0,
            "Maximum temperature should be greater than minimum temperature."
        );
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
    /// higher.
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

    /// Whether to count the calculation of non-cached ICC values as a "step", rather than
    /// the calculation of the whole MCM.
    ///
    /// Enabling this option could lead to more efficient resource utilization.
    pub fn set_count_non_cached(mut self, count_non_cached: bool) -> Self {
        self.count_non_cached = count_non_cached;
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

    /// Sets whether MCMs can attempt to swap between arbitrary pools, or only
    /// their direct neighbors.
    pub fn set_arbitrary_swaps(mut self, arbitrary: bool) -> Self {
        self.arbitrary_swaps = arbitrary;
        self
    }

    /// Attach a sender for a `Collector` to this solver.
    pub fn set_sender(mut self, sender: Sender<CollectorEvent<SolverEvent>>) -> Self {
        self.sender = Some(sender);
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
            *current = current.mutate(rng).0;
            let new_log_e = current.par_log_e(&self.dataset, log_e_cache);

            if new_log_e.total_cmp(&old_log_e).is_lt() {
                deltas_log_e_regressions.push(new_log_e.value - old_log_e.value);
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
        (0..self.pool_amount)
            .map(|i| {
                self.min_temperature
                    + (max_temp - self.min_temperature) * (i as f64)
                        / ((self.pool_amount - 1) as f64)
            })
            .collect()
    }

    fn inverse_linear_temperatures(&self, max_temp: f64) -> Vec<f64> {
        (0..self.pool_amount)
            .map(|i| {
                1.0 / max_temp
                    + (1.0 / self.min_temperature - 1.0 / max_temp) * (i as f64)
                        / ((self.pool_amount - 1) as f64)
            })
            .map(|t| 1.0 / t)
            .collect()
    }

    fn geometric_temperatures(&self, max_temp: f64) -> Vec<f64> {
        let r = (max_temp / self.min_temperature).powf(1.0 / ((self.pool_amount - 1) as f64));

        (0..(self.pool_amount as i32))
            .map(|i| self.min_temperature * r.powi(i))
            .collect()
    }

    fn exponential_temperatures(&self, max_temp: f64) -> Vec<f64> {
        (0..(self.pool_amount as i32))
            .map(|i| {
                (self.min_temperature.ln()
                    + (max_temp.ln() - self.min_temperature.ln()) * (i as f64)
                        / ((self.pool_amount - 1) as f64))
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

    fn get_swaps(&self, turn: usize) -> Vec<(usize, usize)> {
        let mut output = vec![];
        if self.arbitrary_swaps {
            let mut rng = rand::rng();
            let mut pools: Vec<_> = (0..self.pool_amount).collect();
            pools.shuffle(&mut rng);
            output = pools.chunks_exact(2).map(|a| (a[0], a[1])).collect();
        } else {
            if turn.is_multiple_of(2) {
                for i in 0..self.pool_amount / 2 {
                    output.push((i * 2, i * 2 + 1));
                }
            } else {
                for i in 0..self.pool_amount / 2 {
                    if i * 2 + 2 >= self.pool_amount {
                        break;
                    };
                    output.push((i * 2 + 1, i * 2 + 2));
                }
            }
        }
        output
    }
}

impl Solver for ParallelTemperingSolver {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(ParallelTemperingSolver {
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
            arbitrary_swaps: false,
            count_non_cached: false,
            sender: None,
        })
    }

    fn solve(&self) -> SolverReport {
        // let _par_span = info_span!("parallel_tempering").entered();
        // info!(
        //     target: "init",
        //     amount = self.pool_amount,
        //     "Started Parallel Tempering Solver.",
        // );
        let log_e_cache = get_par_log_e_cache();

        let mut starter = match self.starter {
            AnnealingStarter::Single => {
                MinimallyComplexModel::full(NonZero::new(self.dataset.variables()).unwrap())
            }
            AnnealingStarter::Trivial => {
                MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap())
            }
        };

        // let temp_span = info_span!("Initialize Temperature").entered();
        let temperatures = self.calculate_temperatures(&log_e_cache, &mut starter);
        if let Some(tx) = &self.sender {
            tx.send(ParTempData::InitPoolTemps(temperatures.clone()).into())
                .unwrap();
        }

        let mut pools: Vec<_> = temperatures
            .iter()
            .enumerate()
            .map(|(nr, t)| {
                TemperPool::new(
                    starter.clone(),
                    *t,
                    &self.dataset,
                    nr,
                    self.count_non_cached,
                    self.sender.clone(),
                )
            })
            .collect();

        let bar = Arc::new(Mutex::new(par_tqdm!(
            total = self.shuffles * self.steps_per_shuffle * temperatures.len()
        )));

        for i in 0..self.shuffles {
            if let Some(tx) = &self.sender {
                tx.send(ParTempData::NextIteration(i).into()).unwrap();
                let ids = pools.iter().map(|p| p.mcm_id).collect();
                tx.send(ParTempData::MCMPoolIds(ids).into()).unwrap()
            }

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

            let swaps = self.get_swaps(i);
            let rng = rand::rng();
            swap_pools(&mut pools, swaps, rng);
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

fn swap_pools(pools: &mut [TemperPool], swaps: Vec<(usize, usize)>, mut rng: ThreadRng) {
    for (i, j) in swaps {
        let highest = i.max(j);
        let lowest = i.min(j);
        let (start, end) = pools.split_at_mut(highest);
        let p_one = start.get_mut(lowest).unwrap();
        let p_two = end.get_mut(0).unwrap();
        p_one.swap(p_two, &mut rng);
    }
}

pub struct TemperPool {
    mcm: MinimallyComplexModel,
    log_e: f64,
    temperature: f64,
    best_mcm: MinimallyComplexModel,
    best_log_e: f64,
    id: usize,
    mcm_id: usize,
    count_new_iccs: bool,
    sender: Option<Sender<CollectorEvent<SolverEvent>>>,
}

impl TemperPool {
    pub fn new(
        mcm: MinimallyComplexModel,
        temperature: f64,
        dataset: &VecDataset,
        id: usize,
        count_new_iccs: bool,
        sender: Option<Sender<CollectorEvent<SolverEvent>>>,
    ) -> Self {
        let log_e = mcm.log_e(dataset, &mut None);
        Self {
            best_mcm: mcm.clone(),
            best_log_e: log_e,
            mcm,
            log_e,
            temperature,
            id,
            mcm_id: id,
            count_new_iccs,
            sender,
        }
    }

    /// Perform a Metropolis step on the held `MinimallyComplexModel`.
    fn step(
        &mut self,
        dataset: &VecDataset,
        rng: &mut ThreadRng,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
    ) -> usize {
        let (new_mcm, mut_type) = self.mcm.mutate(rng);
        let result = new_mcm.par_log_e(dataset, log_e_cache);
        let new_log_e = result.value;

        let accepted = new_log_e.total_cmp(&self.log_e).is_gt()
            || rng.random_bool(((new_log_e - self.log_e) / self.temperature).exp());
        if accepted {
            if new_log_e.total_cmp(&self.best_log_e).is_gt() {
                self.best_log_e = new_log_e;
                self.best_mcm = new_mcm.clone();
            }

            self.mcm = new_mcm;
            self.log_e = new_log_e;
        }

        if let Some(tx) = &self.sender {
            tx.send(MCMData::Mutation(MutationEvent { mut_type, accepted }).into())
                .unwrap();
        }

        result.new_icc_count
    }

    /// Perform n Metropolis steps. Also updates the progress bar.
    fn step_n(
        &mut self,
        dataset: &VecDataset,
        n: usize,
        log_e_cache: &Option<Arc<DashMap<FixedBitSet, f64>>>,
        bar: &Arc<Mutex<Bar>>,
    ) {
        let mut steps = 0;
        while steps < n {
            let mut rng = rand::rng();
            let new_iccs = self.step(dataset, &mut rng, log_e_cache);
            if self.count_new_iccs {
                steps += new_iccs;
                bar.lock().unwrap().update(new_iccs).unwrap();
            } else {
                steps += 1;
                bar.lock().unwrap().update(1).unwrap();
            }
        }
    }

    /// Stochastically swaps the `MinimallyComplesModel`s between this and the supplied
    /// `TemperPool`.
    fn swap(&mut self, other: &mut TemperPool, rng: &mut ThreadRng) {
        let probability = ((1.0 / other.temperature - 1.0 / self.temperature)
            * (other.log_e - self.log_e))
            .exp()
            .min(1.0);
        let accepted = rng.random_bool(probability);
        if accepted {
            mem::swap(&mut self.mcm, &mut other.mcm);
            mem::swap(&mut self.log_e, &mut other.log_e);
            mem::swap(&mut self.mcm_id, &mut other.mcm_id);
        }
        if let Some(tx) = &self.sender {
            let swap_event = SwapEvent {
                first_pool: self.id as u64,
                second_pool: other.id as u64,
                probability,
                accepted,
            };
            tx.send(ParTempData::Swap(swap_event).into()).unwrap();
        }
    }
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SwapEvent {
    first_pool: u64,
    second_pool: u64,
    probability: f64,
    accepted: bool,
}

#[derive(Debug, Clone)]
pub enum ParTempData {
    InitPoolTemps(Vec<f64>),
    Swap(SwapEvent),
    MCMPoolIds(Vec<usize>),
    NextIteration(usize),
}

impl From<ParTempData> for CollectorEvent<SolverEvent> {
    fn from(value: ParTempData) -> Self {
        CollectorEvent::Data(SolverEvent::ParTemp(value))
    }
}
