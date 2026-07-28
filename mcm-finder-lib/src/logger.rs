use std::{
    default,
    fs::File,
    ops::{Deref, DerefMut},
    process::id,
    sync::{Arc, Mutex, RwLock, mpsc::Sender},
};

use annolog::*;
use csv::Writer;
use num_enum::FromPrimitive;

use crate::solvers::{ParTempData, SwapEvent};

#[derive(Debug, Clone)]
pub enum SolverEvent {
    ParTemp(ParTempData),
    SomethingElse,
}

struct ParTempHandler {
    iteration: usize,
    swap_history_file: Writer<File>,
    mcm_id_file: Writer<File>,
}

impl Default for ParTempHandler {
    fn default() -> Self {
        ParTempHandler {
            iteration: 0,
            swap_history_file: csv::Writer::from_path("./results/swap_history.csv").unwrap(),
            mcm_id_file: csv::Writer::from_path("./results/mcm_ids.csv").unwrap(),
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
            }
        }
    }
}

pub fn get_collector() -> (Collector<SolverEvent>, Sender<CollectorEvent<SolverEvent>>) {
    CollectorBuilder::<SolverEvent>::new()
        .with_handler(ParTempHandler::default())
        .build()
}
