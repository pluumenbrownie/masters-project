use std::{collections::HashMap, num::NonZero, path::Path};

use kdam::{BarExt, tqdm};
use rand::RngExt;

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::{
        anneal_temps::AnnealingTemperature,
        get_log_e_cache,
        solvers_base::{Solver, SolverReport},
    },
};

pub mod anneal_temps;

#[derive(Debug, Default)]
pub enum AnnealingStarter {
    Single,
    #[default]
    Trivial,
}

pub struct SimulatedAnnealingSearcher {
    dataset: VecDataset,
    starter: AnnealingStarter,
    temperature: AnnealingTemperature,
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
        let mut log_e_cache = get_log_e_cache();

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
