use std::{
    cell::{LazyCell, OnceCell},
    default,
    fs::File,
    ops::{Deref, DerefMut},
    process::id,
    sync::{Arc, Mutex, RwLock, mpsc::Sender},
};

use annolog::*;
use csv::Writer;
use num_enum::FromPrimitive;

use crate::{
    mcm::{MCMData, MutationEvent},
    solvers::{AnnealingData, GreedyData, ParTempData, SwapEvent},
};

#[derive(Debug, Clone)]
pub enum SolverEvent {
    ParTemp(ParTempData),
    Greedy(GreedyData),
    Mcm(MCMData),
    Annealing(AnnealingData),
    SomethingElse,
}

pub type SolverEventSender = Sender<CollectorEvent<SolverEvent>>;
// pub type Sendable = impl Into<CollectorEvent<SolverEvent>>;

struct ParTempHandler {
    iteration: usize,
    swap_history_file: LazyCell<Writer<File>>,
    mcm_id_file: LazyCell<Writer<File>>,
    mcm_log_e_file: LazyCell<Writer<File>>,
}

impl Default for ParTempHandler {
    fn default() -> Self {
        ParTempHandler {
            iteration: 0,
            swap_history_file: LazyCell::new(|| {
                csv::Writer::from_path("./results/swap_history.csv").unwrap()
            }),
            mcm_id_file: LazyCell::new(|| csv::Writer::from_path("./results/mcm_ids.csv").unwrap()),
            mcm_log_e_file: LazyCell::new(|| {
                csv::Writer::from_path("./results/par_temp_log_e.csv").unwrap()
            }),
        }
    }
}

impl Handler<SolverEvent> for ParTempHandler {
    fn handle(&mut self, event: &SolverEvent) {
        if let SolverEvent::ParTemp(par_temp_data) = event {
            match par_temp_data {
                ParTempData::InitPoolTemps(temps) => {
                    self.mcm_id_file
                        .write_record(temps.iter().map(|i| format!("{i}")))
                        .unwrap();
                }
                ParTempData::Swap(swap) => {
                    self.swap_history_file.serialize(swap).unwrap();
                    self.swap_history_file.flush().unwrap();
                }
                ParTempData::MCMPoolIds(ids) => {
                    self.mcm_id_file
                        .write_record(ids.iter().map(|i| format!("{i}")))
                        .unwrap();
                    self.mcm_id_file.flush().unwrap();
                }
                ParTempData::NextIteration(iteration) => self.iteration = *iteration,
                ParTempData::LogE(log_e_with_id) => {
                    self.mcm_log_e_file.serialize(log_e_with_id).unwrap();
                    self.mcm_log_e_file.flush().unwrap();
                }
            }
        }
    }
}

struct McmHandler {
    iteration: usize,
    mutation_history_file: LazyCell<Writer<File>>,
    log_e_file: LazyCell<Writer<File>>,
}

impl Default for McmHandler {
    fn default() -> Self {
        McmHandler {
            iteration: 0,
            mutation_history_file: LazyCell::new(|| {
                csv::Writer::from_path("./results/mutation_history.csv").unwrap()
            }),
            log_e_file: LazyCell::new(|| {
                let mut writer = csv::Writer::from_path("./results/log_e_history.csv").unwrap();
                let header = ["log_e"];
                writer.write_record(header).unwrap();
                writer
            }),
        }
    }
}

impl Handler<SolverEvent> for McmHandler {
    fn handle(&mut self, event: &SolverEvent) {
        if let SolverEvent::Mcm(mcm_data) = event {
            match mcm_data {
                MCMData::Mutation(mut_event) => {
                    self.mutation_history_file.serialize(mut_event).unwrap();
                    self.mutation_history_file.flush().unwrap();
                }
                MCMData::LogE(log_e) => {
                    self.log_e_file.write_record([format!("{log_e}")]).unwrap();
                    self.log_e_file.flush().unwrap();
                }
            }
        }
    }
}

struct GreedyHandler {
    iteration: usize,
    log_e_file: LazyCell<Writer<File>>,
}

impl Default for GreedyHandler {
    fn default() -> Self {
        GreedyHandler {
            iteration: 0,
            log_e_file: LazyCell::new(|| {
                csv::Writer::from_path("./results/log_e_greedy.csv").unwrap()
            }),
        }
    }
}

impl Handler<SolverEvent> for GreedyHandler {
    fn handle(&mut self, event: &SolverEvent) {
        if let SolverEvent::Greedy(greedy_data) = event {
            match greedy_data {
                GreedyData::LogE(greedy_log_e) => {
                    self.log_e_file.serialize(greedy_log_e).unwrap();
                    self.log_e_file.flush().unwrap();
                }
            }
        }
    }
}

struct AnnealingHandler {
    iteration: usize,
    log_e_file: LazyCell<Writer<File>>,
}

impl Default for AnnealingHandler {
    fn default() -> Self {
        AnnealingHandler {
            iteration: 0,
            log_e_file: LazyCell::new(|| {
                csv::Writer::from_path("./results/log_e_annealing.csv").unwrap()
            }),
        }
    }
}

impl Handler<SolverEvent> for AnnealingHandler {
    fn handle(&mut self, event: &SolverEvent) {
        if let SolverEvent::Annealing(annealing_data) = event {
            match annealing_data {
                AnnealingData::LogE(annealing_log_e) => {
                    self.log_e_file.serialize(annealing_log_e).unwrap();
                    self.log_e_file.flush().unwrap();
                }
            }
        }
    }
}

pub fn get_collector() -> (Collector<SolverEvent>, SolverEventSender) {
    CollectorBuilder::<SolverEvent>::new()
        .with_handler(ParTempHandler::default())
        .with_handler(McmHandler::default())
        .with_handler(GreedyHandler::default())
        .with_handler(AnnealingHandler::default())
        .build()
}
