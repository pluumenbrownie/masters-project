use kdam::{BarExt, tqdm};
use mcm_finder_lib::{dataset::Dataset, mcm::MinimallyComplexModel};
use std::{collections::HashMap, num::NonZero, path::Path};

fn main() {
    // let dataset = Dataset::read_from_file(Path::new(
    //     "mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat",
    // ))
    // .unwrap();
    let dataset =
        Dataset::read_from_file(Path::new("mcm-finder-lib/tests/data/MNIST11.sorted")).unwrap();
    // let dataset =
    //     Dataset::read_from_file(Path::new("mcm-finder-lib/tests/data/Big5PT.sorted")).unwrap();
    let variables = dataset.variables();
    println!("Greedy Searching for MCM of dataset with");
    println!("  - {variables} variables");
    println!("  - {} datapoints", dataset.datapoints());
    println!("  - {} bins", dataset.bins());

    let mut best_mcm = MinimallyComplexModel::trivial(NonZero::new(variables).unwrap());

    // let mut best_mcm: MinimallyComplexModel = current.clone();
    let mut log_e_cache = Some(HashMap::new());
    let mut best_log_e = best_mcm.log_e(&dataset, &mut log_e_cache);

    let mut length = 0usize;
    for parts_left in (0usize..variables).rev() {
        for basis in tqdm!(1..=parts_left) {
            for _ in 0..basis {
                length += 1;
            }
        }
    }

    // we merge one partition each round
    let mut progress = tqdm!(total = length);
    for parts_left in (0usize..variables).rev() {
        let mut new_best = false;
        let old_best = best_mcm.clone();
        progress.set_description(format!("{parts_left} steps left"));
        for basis in 1..=parts_left {
            for into in 0..basis {
                // generation.push(best_mcm.merge(basis, into));
                let candidate = old_best.merge(basis, into);
                if candidate
                    .log_e(&dataset, &mut log_e_cache)
                    .total_cmp(&best_log_e)
                    .is_ge()
                {
                    best_mcm = candidate.clone();
                    best_log_e = best_mcm.log_e(&dataset, &mut log_e_cache);
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

    let mcm = best_mcm;
    println!("Best MCM with log E {}", best_log_e);
    println!("{}", mcm);
}
