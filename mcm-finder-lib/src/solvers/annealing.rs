use std::{collections::HashMap, num::NonZero, path::Path};

use annolog::CollectorEvent;
use kdam::{BarExt, tqdm};
use rand::RngExt;
use serde::Serialize;

use crate::{
    dataset::{Dataset, simple::VecDataset},
    logger::{SolverEvent, SolverEventSender},
    mcm::{MCMData, MinimallyComplexModel, MutationEvent},
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

pub struct SimulatedAnnealingSolver<T> {
    dataset: T,
    starter: AnnealingStarter,
    temperature: AnnealingTemperature,
    sender: Option<SolverEventSender>,
}

impl<T: Dataset> SimulatedAnnealingSolver<T> {
    pub fn set_starter(mut self, starter: AnnealingStarter) -> Self {
        self.starter = starter;
        self
    }

    pub fn set_temperature(mut self, temperature: AnnealingTemperature) -> Self {
        self.temperature = temperature;
        self
    }

    /// Attach a sender for a `Collector` to this solver.
    pub fn set_sender(mut self, sender: SolverEventSender) -> Self {
        self.sender = Some(sender);
        self
    }
}

impl<T: Dataset> Solver for SimulatedAnnealingSolver<T> {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(SimulatedAnnealingSolver {
            dataset: T::read_from_file(filepath)?,
            starter: AnnealingStarter::default(),
            temperature: AnnealingTemperature::default(),
            sender: None,
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
            let (candidate, mut_type) = current.mutate(&mut rng);
            let candidate_log_e = candidate.log_e(&self.dataset, &mut log_e_cache);
            let current_log_e = current.log_e(&self.dataset, &mut log_e_cache);

            if candidate_log_e.total_cmp(&best_log_e).is_gt() {
                best_mcm = candidate.clone();
                best_log_e = candidate.log_e(&self.dataset, &mut log_e_cache);
            }
            let mut accepted = false;
            if candidate_log_e.total_cmp(&current_log_e).is_gt()
                || rng.random_bool(((candidate_log_e - current_log_e) / temp).exp())
            {
                current = candidate;
                accepted = true;
            }

            progress.set_description(format!("T={:.3}/{} - {:.1}", temp, end, best_log_e));
            self.send(MCMData::Mutation(MutationEvent {
                mut_type,
                accepted,
                temperature: temp,
            }))
            .unwrap();
            self.send(AnnealingData::LogE(AnnealingLogE {
                log_e: current.log_e(&self.dataset, &mut log_e_cache),
                temperature: temp,
            }))
            .unwrap();
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

    fn get_sender(&self) -> Option<&SolverEventSender> {
        self.sender.as_ref()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct AnnealingLogE {
    pub log_e: f64,
    pub temperature: f64,
}

#[derive(Debug, Clone)]
pub enum AnnealingData {
    LogE(AnnealingLogE),
}

impl From<AnnealingData> for CollectorEvent<SolverEvent> {
    fn from(value: AnnealingData) -> Self {
        CollectorEvent::Data(SolverEvent::Annealing(value))
    }
}
