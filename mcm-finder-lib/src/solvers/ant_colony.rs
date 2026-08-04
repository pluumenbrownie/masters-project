use std::{
    collections::HashMap,
    env::var,
    ffi::FromBytesUntilNulError,
    fmt::{Debug, Display},
    ops::Deref,
    sync::{Arc, Mutex},
};

use kdam::{BarExt, par_tqdm};
use rand::{
    rngs::ThreadRng,
    seq::{IndexedRandom, index::sample_weighted},
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::distribution::ChiError;

use crate::{
    dataset::Dataset,
    mcm::{self, MinimallyComplexModel},
    solvers::{ConstructiveStrategy::VarianceFirst, Solver, SolverReport, get_par_log_e_cache},
};

#[derive(Debug, Clone, Copy)]
enum ICCorDivider {
    Icc(usize),
    Divider(usize),
}

#[derive(Debug, Clone)]
struct AntPath {
    data: Vec<ICCorDivider>,
    variables: usize,
}

impl From<&AntPath> for MinimallyComplexModel {
    fn from(permutation: &AntPath) -> Self {
        let mut mcm = vec![0usize; permutation.variables];

        let mut icc_counter = 1usize;
        for x in &permutation.data {
            match x {
                ICCorDivider::Divider(_) => icc_counter += 1,
                ICCorDivider::Icc(variable) => mcm[*variable] = icc_counter,
            }
        }

        MinimallyComplexModel::from_vector(mcm)
    }
}

impl AntPath {
    fn generate(pheromone_environment: &PheromoneEnvironment, rng: &mut ThreadRng) -> AntPath {
        let mut permutation = vec![ICCorDivider::Icc(0)];
        let mut last_entry = ICCorDivider::Icc(0);
        let mut destinations: Vec<ICCorDivider> = (1..pheromone_environment.nodes)
            .map(|i| {
                if i > pheromone_environment.variables() - 1 {
                    ICCorDivider::Divider(i - pheromone_environment.variables())
                } else {
                    ICCorDivider::Icc(i)
                }
            })
            .collect();

        while destinations
            .iter()
            .any(|i| matches!(i, ICCorDivider::Icc(_)))
        {
            let get_weight = |to| pheromone_environment.get(last_entry, destinations[to]);
            let choice_index = sample_weighted(rng, destinations.len(), get_weight, 1)
                .unwrap()
                .index(0);
            let choice = destinations[choice_index];
            if let ICCorDivider::Icc(_) = choice {
                destinations.swap_remove(choice_index);
            }
            last_entry = choice;
            permutation.push(choice);
        }

        AntPath {
            data: permutation,
            variables: pheromone_environment.variables(),
        }
    }
}

struct PheromoneEnvironment {
    data: Vec<f64>,
    nodes: usize,
    evaporation_rate: f64,
    minimum: f64,
    excretion_amount: f64,
}

impl PheromoneEnvironment {
    fn new(variables: usize) -> PheromoneEnvironment {
        // 0..variables: mcm variables; variables..2*variables: divider;
        let nodes = variables * 2;
        let mut data = vec![0.0f64; nodes.pow(2)];

        for (nr, p) in data.iter_mut().enumerate() {
            let from = nr % nodes;
            let to = nr / nodes;

            if from != to && (from < variables && to < variables) {
                *p = 1.0;
            } else if from != to && (from < variables || to < variables) {
                *p = 0.01;
            }
        }

        PheromoneEnvironment {
            data,
            nodes,
            evaporation_rate: 0.0,
            minimum: 0.0,
            excretion_amount: 0.0,
        }
    }

    fn variables(&self) -> usize {
        self.nodes / 2
    }

    fn set_params(mut self, evaporation_rate: f64, minimum: f64, excretion_amount: f64) -> Self {
        self.evaporation_rate = evaporation_rate;
        self.minimum = minimum;
        self.excretion_amount = excretion_amount;
        self
    }

    fn decode_enum(&self, from: ICCorDivider) -> usize {
        match from {
            ICCorDivider::Icc(value) => value,
            ICCorDivider::Divider(value) => self.nodes / 2 + value,
        }
    }

    fn get(&self, from: ICCorDivider, to: ICCorDivider) -> f64 {
        let from = self.decode_enum(from);
        let to = self.decode_enum(to);
        self.data[from + to * (self.nodes)]
    }

    fn get_mut(&mut self, from: ICCorDivider, to: ICCorDivider) -> &mut f64 {
        let from = self.decode_enum(from);
        let to = self.decode_enum(to);
        self.data
            .get_mut(from + to * (self.nodes))
            .expect("Bad ends received.")
    }

    fn get_weights(&self, from: ICCorDivider, to_slice: &[ICCorDivider]) -> Vec<f64> {
        let from = self.decode_enum(from);

        to_slice
            .iter()
            .map(|to| {
                let to = self.decode_enum(*to);
                self.data[from + to * (self.nodes)]
            })
            .collect()
    }

    fn evaporate(&mut self) {
        for x in self.data.iter_mut() {
            if *x != 0.0 {
                *x = self.minimum.max(*x - *x * self.evaporation_rate);
            }
        }
    }

    fn excrete(&mut self, path: AntPath) {
        for ends in path.data.windows(2) {
            let from = ends[0];
            let to = ends[1];
            *self.get_mut(from, to) += self.excretion_amount;
        }
    }
}

impl Debug for PheromoneEnvironment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("PheromoneEnvironment")
            .field(&self.nodes)
            .field(&self.evaporation_rate)
            .field(&self.minimum)
            .field(&self.excretion_amount)
            .finish_non_exhaustive()?;
        write!(f, "\ndata: ")?;
        for (nr, x) in self.data.iter().enumerate() {
            if nr % self.nodes == 0 {
                writeln!(f)?;
            }
            write!(f, "{x:.2} ")?;
        }
        Ok(())
    }
}

pub struct AntColonyOptimization<D: Dataset> {
    dataset: D,
    steps: usize,
    ants: usize,
    amount_update: usize,
    evaporation_rate: f64,
    minimum: f64,
    excretion_amount: f64,
}

impl<D: Dataset + Sync> Solver for AntColonyOptimization<D> {
    fn from_file(filepath: &std::path::Path) -> Result<Self, crate::mcm_error::MCMError>
    where
        Self: Sized,
    {
        Ok(AntColonyOptimization {
            dataset: D::read_from_file(filepath)?,
            steps: 100,
            ants: 160,
            amount_update: 30,
            evaporation_rate: 0.1,
            minimum: 0.01,
            excretion_amount: 0.02,
        })
    }

    fn solve(&self) -> super::SolverReport {
        let log_e_cache = get_par_log_e_cache();
        let mut mut_rng = rand::rng();
        let rng = &mut mut_rng;
        let mut best_mcm_option = None;

        let bar = Arc::new(Mutex::new(par_tqdm!(total = self.ants * self.steps)));

        let mut environment = PheromoneEnvironment::new(self.dataset.variables()).set_params(
            self.evaporation_rate,
            self.minimum,
            self.excretion_amount,
        );

        for _ in 0..self.steps {
            let ant_paths: Vec<AntPath> = (0..self.ants)
                .map(|_| AntPath::generate(&environment, rng))
                .collect();
            let mut log_e_pairs: Vec<_> = ant_paths
                .into_par_iter()
                .map(|path| {
                    let mcm = MinimallyComplexModel::from(&path);
                    let log_e = mcm.par_log_e(&self.dataset, &log_e_cache);
                    bar.lock().unwrap().update(1).unwrap();
                    (path, log_e)
                })
                .collect();

            log_e_pairs.sort_by(|a, b| a.1.total_cmp(&b.1).reverse());
            best_mcm_option = Some(log_e_pairs[0].clone());

            environment.evaporate();
            for (path, _) in log_e_pairs.into_iter().take(self.amount_update) {
                environment.excrete(path);
            }
        }

        println!("{:?}", environment);

        let best_mcm = best_mcm_option.unwrap();
        SolverReport::new(
            MinimallyComplexModel::from(&best_mcm.0),
            *best_mcm.1,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{:.0}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}
