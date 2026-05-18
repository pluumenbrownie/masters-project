use std::{cell::RefCell, collections::HashMap, num::NonZero, path::Path};

use fixedbitset::FixedBitSet;
use kdam::{BarExt, tqdm};
use rand::{RngExt, rngs::ThreadRng};

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::{
        AnnealingStarter, get_log_e_cache,
        solvers_base::{Solver, SolverReport},
    },
};

pub struct AdaptiveAnnealingSearcher {
    dataset: VecDataset,
    starter: AnnealingStarter,
    temperature: RefCell<AdaptiveTemperature>,
}

impl AdaptiveAnnealingSearcher {
    pub fn set_starter(mut self, starter: AnnealingStarter) -> Self {
        self.starter = starter;
        self
    }

    fn calculate_inital_temperature(
        &self,
        current: &mut MinimallyComplexModel,
        rng: &mut ThreadRng,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
    ) {
        let mut deltas_log_e_regressions: Vec<f64> = vec![];

        while deltas_log_e_regressions.len() < 100 {
            let old_log_e = current.log_e(&self.dataset, log_e_cache);
            *current = current.mutate(rng);
            let new_log_e = current.log_e(&self.dataset, log_e_cache);

            if new_log_e.total_cmp(&old_log_e).is_lt() {
                deltas_log_e_regressions.push(new_log_e - old_log_e);
            }
        }

        self.temperature
            .borrow_mut()
            .calculate_start_temperature(deltas_log_e_regressions);
    }
}

impl Solver for AdaptiveAnnealingSearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(AdaptiveAnnealingSearcher {
            dataset: VecDataset::read_from_file(filepath)?,
            starter: AnnealingStarter::default(),
            temperature: AdaptiveTemperature::new(0.1, 100).into(),
        })
    }

    fn solve(&self) -> SolverReport {
        let mut current = match self.starter {
            AnnealingStarter::Single => {
                MinimallyComplexModel::full(NonZero::new(self.dataset.variables()).unwrap())
            }
            AnnealingStarter::Trivial => {
                MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap())
            }
        };

        let mut rng = rand::rng();
        let mut log_e_cache = get_log_e_cache();
        let mut steps_since_current_improvement = 0usize;

        let mut best_mcm = current.clone();
        let mut best_log_e = current.log_e(&self.dataset, &mut log_e_cache);

        let mut progress = tqdm!();
        self.calculate_inital_temperature(&mut current, &mut rng, &mut log_e_cache);

        // while temp > self.temperature.end {
        let mut adaptive_temperature_iter = self.temperature.borrow_mut().into_iter();

        while let Some(temp) = adaptive_temperature_iter.next() {
            let candidate = current.mutate(&mut rng);
            let candidate_log_e = candidate.log_e(&self.dataset, &mut log_e_cache);
            let current_log_e = current.log_e(&self.dataset, &mut log_e_cache);

            if candidate_log_e.total_cmp(&best_log_e).is_gt() {
                best_mcm = candidate.clone();
                best_log_e = candidate.log_e(&self.dataset, &mut log_e_cache);
            }
            if candidate_log_e.total_cmp(&current_log_e).is_gt()
                || rng.random_bool(((candidate_log_e - current_log_e) / temp).exp())
            {
                current = candidate;
                steps_since_current_improvement = 0;
            } else {
                steps_since_current_improvement += 1;
            }

            if steps_since_current_improvement > 100 {
                steps_since_current_improvement = 0;
                adaptive_temperature_iter.reset();
            }

            progress.set_description(format!(
                "T={:.3} | Best Log(E)={:.0} | Current Log(E)={:.0}",
                temp, best_log_e, current_log_e
            ));
            let _ = progress.update(1);
        }
        SolverReport::new(
            best_mcm,
            best_log_e,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}

#[derive(Default, Clone, Copy)]
pub struct AdaptiveTemperature {
    start: f64,
    step: usize,
    temp: f64,
    end: f64,
    const_steps: usize,
}

impl AdaptiveTemperature {
    pub fn new(end: f64, const_steps: usize) -> Self {
        Self {
            end,
            const_steps,
            ..Default::default()
        }
    }

    fn calculate_start_temperature(&mut self, deltas_log_e_regressions: Vec<f64>) {
        let average_regression =
            deltas_log_e_regressions.iter().sum::<f64>() / deltas_log_e_regressions.len() as f64;
        self.start = average_regression / 0.9_f64.ln();
    }

    fn delta(&mut self) {
        self.temp = self.start
            * 0.997f64.powi(
                (self.step.saturating_sub(self.const_steps))
                    .try_into()
                    .unwrap(),
            );
    }

    fn reset(&mut self) {
        self.start *= 0.8;
        self.temp = self.start;
        self.step = 0;
    }
}

impl Iterator for AdaptiveTemperature {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        self.delta();
        self.step += 1;
        let output = Some(self.temp);
        if self.start < self.end || self.temp < f64::EPSILON {
            None
        } else {
            output
        }
    }
}
