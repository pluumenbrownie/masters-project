use std::{collections::HashMap, marker::Sized, path::Path};

use fixedbitset::FixedBitSet;
use kdam::tqdm;

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    mcm_error::MCMError,
    solvers::solvers_base::{Solver, SolverReport},
};

pub struct ExhaustiveSearcher {
    pub(crate) dataset: VecDataset,
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

pub(crate) fn generate_all_mcms(number: usize) -> Vec<MinimallyComplexModel> {
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

pub(crate) fn add_to(current: &mut Vec<usize>, index: usize, output: &mut Vec<Vec<usize>>) {
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
