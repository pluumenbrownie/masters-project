use std::{
    collections::{HashMap, VecDeque},
    num::NonZero,
    path::Path,
};

use annolog::CollectorEvent;
use kdam::{BarExt, tqdm};
use rand::RngExt;
use serde::Serialize;

use crate::{
    dataset::{Dataset, simple::VecDataset},
    logger::{SolverEvent, SolverEventSender},
    mcm::{MCMData::LogE, MinimallyComplexModel},
    mcm_error::MCMError,
    solvers::{
        AnnealingStarter, get_log_e_cache,
        solvers_base::{Solver, SolverReport},
    },
};

pub struct HillClimberSolver<T> {
    dataset: T,
    starter: AnnealingStarter,
    max_steps: usize,
    history_size: usize,
    stagnation_steps: usize,
    sender: Option<SolverEventSender>,
}

impl<T: Dataset> HillClimberSolver<T> {
    pub fn set_starter(mut self, starter: AnnealingStarter) -> Self {
        self.starter = starter;
        self
    }

    pub fn set_max_steps(mut self, max_steps: usize) -> Self {
        assert!(max_steps > 0, "max_steps must be greater than zero.");
        self.max_steps = max_steps;
        self
    }

    pub fn set_history_size(mut self, size: usize) -> Self {
        assert!(size > 0, "max_steps must be greater than zero.");
        self.history_size = size;
        self
    }

    pub fn set_stagnation_steps(mut self, steps: usize) -> Self {
        self.stagnation_steps = steps;
        self
    }

    /// Attach a sender for a `Collector` to this solver.
    pub fn set_sender(mut self, sender: SolverEventSender) -> Self {
        self.sender = Some(sender);
        self
    }
}

impl<T: Dataset> Solver for HillClimberSolver<T> {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(HillClimberSolver {
            dataset: T::read_from_file(filepath)?,
            starter: AnnealingStarter::default(),
            max_steps: 10_000,
            history_size: 1,
            stagnation_steps: 1000,
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
        let mut history = VecDeque::with_capacity(self.history_size);
        history.push_back(current.clone());
        let mut stagnation_counter = 0usize;

        let mut progress = tqdm!(total = self.max_steps);

        for _ in 0..self.max_steps {
            let (candidate, _) = current.mutate(&mut rng);
            let candidate_log_e = candidate.log_e(&self.dataset, &mut log_e_cache);
            let old_log_e = history
                .front()
                .map(|mcm: &MinimallyComplexModel| mcm.log_e(&self.dataset, &mut log_e_cache));

            if old_log_e.is_some_and(|log_e| candidate_log_e > log_e)
                | (candidate_log_e > current.log_e(&self.dataset, &mut log_e_cache))
            {
                if candidate_log_e > current.log_e(&self.dataset, &mut log_e_cache) {
                    // if candidate_log_e > old_log_e.unwrap() {
                    if history.len() >= self.history_size {
                        history.pop_front();
                    }
                    history.push_back(candidate.clone());
                }

                current = candidate.clone();
                if candidate_log_e > best_log_e {
                    best_log_e = candidate.log_e(&self.dataset, &mut log_e_cache);
                    best_mcm = candidate;
                }
                stagnation_counter = 0;
            } else {
                stagnation_counter += 1;
                if stagnation_counter > self.stagnation_steps {
                    break;
                }
            }

            progress.set_description(format!("Log E: {:.1}", best_log_e));
            self.send(HillClimbingData::LogE(HillClimbingLogE {
                log_e: current.log_e(&self.dataset, &mut log_e_cache),
                historic_log_e: old_log_e.unwrap_or(current.log_e(&self.dataset, &mut log_e_cache)),
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
pub struct HillClimbingLogE {
    pub log_e: f64,
    pub historic_log_e: f64,
}

#[derive(Debug, Clone)]
pub enum HillClimbingData {
    LogE(HillClimbingLogE),
}

impl From<HillClimbingData> for CollectorEvent<SolverEvent> {
    fn from(value: HillClimbingData) -> Self {
        CollectorEvent::Data(SolverEvent::HillClimbing(value))
    }
}
