use std::{collections::HashMap, fmt::Display, marker::Sized, num::NonZero, path::Path};

use fixedbitset::FixedBitSet;
use kdam::{BarExt, tqdm};

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
};
pub mod dataset;
pub mod mcm;
mod mcm_error;

#[derive(Debug, Clone)]
pub struct SolverReport {
    mcm: MinimallyComplexModel,
    log_e: f64,
    other_stuff: HashMap<String, String>,
}

impl SolverReport {
    fn new(
        mcm: MinimallyComplexModel,
        log_e: f64,
        other_stuff: HashMap<String, String>,
    ) -> SolverReport {
        SolverReport {
            mcm,
            log_e,
            other_stuff,
        }
    }
}

impl Display for SolverReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Found MCM with Log E of {} and ICCs:\n{}{}",
            self.log_e,
            self.mcm,
            self.other_stuff
                .iter()
                .flat_map(|(k, v)| format!("\n - {:?}: {:?}", k, v).chars().collect::<Vec<_>>())
                .collect::<String>()
        )
    }
}

pub trait Solver {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized;

    fn solve(&self) -> SolverReport;
}

pub struct ExhaustiveSearcher {
    dataset: VecDataset,
}

impl Solver for ExhaustiveSearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(ExhaustiveSearcher {
            dataset: VecDataset::read_from_file(filepath)?,
        })
    }

    fn solve(&self) -> SolverReport {
        let mcms = generate_all_mcms(self.dataset.variables());
        let mut best_mcm: Option<MinimallyComplexModel> = None;
        let mut best_log_e = f64::NEG_INFINITY;

        for mcm in tqdm!(mcms.into_iter()) {
            let log_e = mcm.log_e(&self.dataset, &mut None);
            if log_e > best_log_e {
                best_log_e = log_e;
                best_mcm = Some(mcm);
            }
        }

        let mcm = best_mcm.unwrap();
        SolverReport::new(mcm, best_log_e, HashMap::new())
    }
}

fn generate_all_mcms(number: usize) -> Vec<MinimallyComplexModel> {
    let mut output: Vec<Vec<usize>> = vec![];
    let mut current = vec![0usize; number];

    add_to(&mut current, 0, &mut output);

    output
        .into_iter()
        .map(|v| {
            let mut partitions = vec![
                FixedBitSet::with_capacity_and_blocks(number, [0b000000000]);
                *v.iter().max().unwrap() + 1
            ];
            for (nr, &i) in v.iter().enumerate() {
                partitions[i].set(nr, true);
            }
            MinimallyComplexModel::from_iccs(partitions).unwrap()
        })
        .collect()
}

fn add_to(current: &mut Vec<usize>, index: usize, output: &mut Vec<Vec<usize>>) {
    loop {
        if index < current.len() - 1 {
            add_to(current, index + 1, output);
        }
        if index == current.len() - 1 {
            output.push(current.clone());
        }
        if index == 0 || current[..index].iter().all(|n| *n < current[index]) {
            break;
        } else {
            current[index] += 1;
            current.splice(index + 1.., vec![0usize; current.len() - index - 1]);
        }
    }
}

pub struct GreedySearcher {
    dataset: VecDataset,
}

impl Solver for GreedySearcher {
    fn from_file(filepath: &Path) -> Result<Self, MCMError>
    where
        Self: Sized,
    {
        Ok(GreedySearcher {
            dataset: VecDataset::read_from_file(filepath)?,
        })
    }

    fn solve(&self) -> SolverReport {
        let mut best_mcm =
            MinimallyComplexModel::trivial(NonZero::new(self.dataset.variables()).unwrap());

        // let mut best_mcm: MinimallyComplexModel = current.clone();
        let mut log_e_cache = Some(HashMap::new());
        let mut best_log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache);

        let mut length = 0usize;
        for parts_left in (0usize..self.dataset.variables()).rev() {
            for basis in tqdm!(1..=parts_left) {
                for _ in 0..basis {
                    length += 1;
                }
            }
        }

        // we merge one partition each round
        let mut progress = tqdm!(total = length);
        for parts_left in (0usize..self.dataset.variables()).rev() {
            let mut new_best = false;
            let old_best = best_mcm.clone();
            progress.set_description(format!("{parts_left} steps left"));
            for basis in 1..=parts_left {
                for into in 0..basis {
                    // generation.push(best_mcm.merge(basis, into));
                    let candidate = old_best.merge(basis, into);
                    if candidate
                        .log_e(&self.dataset, &mut log_e_cache)
                        .total_cmp(&best_log_e)
                        .is_ge()
                    {
                        best_mcm = candidate.clone();
                        best_log_e = best_mcm.log_e(&self.dataset, &mut log_e_cache);
                        new_best = true;
                    }
                    let _ = progress.update(1);
                }
            }
            // let candidate = generation
            //     .iter()
            //     .max_by(|&a, &b| a.log_e(&dataset).total_cmp(&b.log_e(&dataset)))
            //     .unwrap();
            // if new >= old
            if !new_best {
                let _ = progress.update_to(length);
                break;
            }
        }

        SolverReport::new(best_mcm, best_log_e, HashMap::new())
    }
}

// #[derive(Debug, Default)]
// pub struct BasisSet {
//     basis_vectors: Vec<FixedBitSet>,
// }

// impl BasisSet {
//     pub fn new(vectors: Vec<FixedBitSet>) -> BasisSet {
//         BasisSet {
//             basis_vectors: vectors,
//         }
//     }

//     pub fn add(&mut self, basis: FixedBitSet) -> Result<(), MCMError> {
//         if self.basis_vectors.iter().any(|b| !b.is_disjoint(&basis)) {
//             return Err(MCMError::OverlappingBasis);
//         }
//         self.basis_vectors.push(basis);
//         Ok(())
//     }

// pub fn create_kset(&self, dataset: Dataset) -> KSet {
//     let new_vectors: Vec<(FixedBitSet, usize)> = dataset
//         .data
//         .iter()
//         .map(|(i, &n)| self.transform_mu_basis(i, n))
//         .collect();
//     // transformed mu vectors are not guaranteed to be unique, so we want to
//     // add them together without loosing information
//     let mut hashmap = HashMap::new();
//     for (o, n) in new_vectors {
//         hashmap.entry(o).and_modify(|v| *v += n).or_insert(n);
//     }
//     KSet {
//         data: hashmap,
//         datapoints: dataset.datapoints,
//     }
// }

// fn transform_mu_basis(&self, i: &FixedBitSet, n: usize) -> (FixedBitSet, usize) {
//     let mut o = FixedBitSet::with_capacity(self.basis_vectors.len());
//     for (nr, b) in self.basis_vectors.iter().enumerate() {
//         if (i & b).count_ones(..) % 2 == 1 {
//             o.insert(nr);
//         }
//     }
//     (o, n)
// }
// }

// #[derive(Debug)]
// pub struct KSet {
//     pub data: HashMap<FixedBitSet, usize>,
//     pub datapoints: usize,
// }
