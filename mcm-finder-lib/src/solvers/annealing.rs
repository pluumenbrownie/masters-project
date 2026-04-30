use std::{collections::HashMap, num::NonZero, path::Path};

use kdam::{BarExt, tqdm};
use rand::RngExt;

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::solvers_base::{Solver, SolverReport},
};

#[derive(Debug, Default)]
pub enum AnnealingStarter {
    Single,
    #[default]
    Trivial,
}

#[derive(Debug)]
pub struct AnnealingTemperature {
    pub(crate) start: f64,
    pub(crate) end: f64,
    pub(crate) decrease_per_step: f64,
    pub(crate) next_temp: Option<Box<AnnealingTemperature>>,
}

impl Default for AnnealingTemperature {
    fn default() -> Self {
        AnnealingTemperature {
            start: 500.0,
            end: 5.0,
            decrease_per_step: 0.1,
            next_temp: None,
        }
    }
}

impl AnnealingTemperature {
    /// Create a new temperature curve for use in simulated annealing. The temperature
    /// will start at `start` and then be decayed exponentially by removing
    /// `temp * decrease_per_step` each step, until `temp < end`.
    ///
    /// The curve can be extended with an additional curve with different `end` and
    /// `decrease_per_step` values with the `then()` method.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::AnnealingTemperature;
    /// let temp_curve = AnnealingTemperature::new(10_000.0, 5.0, 0.0003).then(0.001, 0.00001);
    /// ```
    pub fn new(start: f64, end: f64, decrease_per_step: f64) -> AnnealingTemperature {
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        assert!(
            decrease_per_step > 0.0 || decrease_per_step < 1.0,
            "Decrease should be number between zero and one."
        );
        AnnealingTemperature {
            start,
            end,
            decrease_per_step,
            next_temp: None,
        }
    }

    pub fn steps(&self) -> usize {
        let steps = (self.end / self.start).log2() / (1.0 - self.decrease_per_step).log2();
        match &self.next_temp {
            Some(t) => steps.ceil() as usize + t.steps(),
            None => steps.ceil() as usize,
        }
    }

    pub(crate) fn set_next_none(&mut self, next_value: AnnealingTemperature) {
        match &mut self.next_temp {
            Some(at) => at.set_next_none(next_value),
            None => self.next_temp = Some(Box::new(next_value)),
        };
    }

    pub(crate) fn get_last_end(&self) -> f64 {
        match &self.next_temp {
            Some(at) => at.get_last_end(),
            None => self.end,
        }
    }

    /// Extend this temperature curve with an additional curve with different `end` and
    /// `decrease_per_step` values.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::AnnealingTemperature;
    /// let temp_curve = AnnealingTemperature::new(10_000.0, 5.0, 0.0003).then(0.001, 0.00001);
    /// ```
    pub fn then(mut self, end: f64, decrease_per_step: f64) -> Self {
        let last_end = self.get_last_end();
        self.set_next_none(AnnealingTemperature::new(last_end, end, decrease_per_step));
        self
    }

    /// Returns iterator over the temperature values for this curve.
    pub fn create_iter(&self) -> AnnealTempIter {
        let mut iterator = AnnealTempIter::new(self.start, self.end, self.decrease_per_step, None);
        iterator.next_temp = self.next_temp.as_ref().map(|t| Box::new(t.create_iter()));
        iterator
    }
}

#[derive(Debug, Clone)]
pub struct AnnealTempIter {
    pub(crate) temp: f64,
    pub(crate) end: f64,
    pub(crate) decrease_per_step: f64,
    pub(crate) next_temp: Option<Box<AnnealTempIter>>,
}

impl AnnealTempIter {
    pub(crate) fn new(
        temp: f64,
        end: f64,
        decrease_per_step: f64,
        next_temp: Option<Box<AnnealTempIter>>,
    ) -> AnnealTempIter {
        AnnealTempIter {
            temp,
            end,
            decrease_per_step,
            next_temp,
        }
    }

    pub fn get_current_target(&self) -> &f64 {
        &self.end
    }
}

impl Iterator for AnnealTempIter {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        let output = if self.temp.total_cmp(&self.end).is_ge() {
            Some((self.temp, self.end))
        } else {
            match &self.next_temp {
                Some(t) => {
                    self.end = t.end;
                    self.decrease_per_step = t.decrease_per_step;
                    self.next_temp = t.next_temp.clone();
                    Some((self.temp, self.end))
                }
                None => None,
            }
        };
        self.temp -= self.temp * self.decrease_per_step;
        output
    }
}

pub struct SimulatedAnnealingSearcher {
    pub(crate) dataset: VecDataset,
    pub(crate) starter: AnnealingStarter,
    pub(crate) temperature: AnnealingTemperature,
}

impl SimulatedAnnealingSearcher {
    pub fn set_starter(mut self, starter: AnnealingStarter) -> Self {
        self.starter = starter;
        self
    }

    pub fn set_temperature(mut self, temperature: AnnealingTemperature) -> Self {
        self.temperature = temperature;
        self
    }
}

impl Solver for SimulatedAnnealingSearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(SimulatedAnnealingSearcher {
            dataset: VecDataset::read_from_file(filepath)?,
            starter: AnnealingStarter::default(),
            temperature: AnnealingTemperature::default(),
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
        let mut log_e_cache = Some(HashMap::new());

        let mut best_mcm = current.clone();
        let mut best_log_e = current.log_e(&self.dataset, &mut log_e_cache);

        // let mut temp = self.temperature.start;
        let steps = self.temperature.steps();
        let mut progress = tqdm!(total = steps);

        // while temp > self.temperature.end {
        for (temp, end) in self.temperature.create_iter() {
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
            }

            progress.set_description(format!("T={:.3}/{} - {:.1}", temp, end, best_log_e));
            let _ = progress.update(1);
        }
        SolverReport::new(best_mcm, best_log_e, HashMap::new())
    }
}
