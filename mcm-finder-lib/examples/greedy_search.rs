use kdam::tqdm;
use mcm_finder_lib::{mcm::MinimallyComplexModel, *};
use std::{num::NonZero, path::Path};

fn main() {
    let dataset = dataset::Dataset::read_from_file(Path::new(
        "mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat",
    ))
    .unwrap();

    let mut best_mcm = MinimallyComplexModel::trivial(NonZero::new(9).unwrap());

    // let mut best_mcm: MinimallyComplexModel = current.clone();
    let mut best_log_e = best_mcm.log_e(&dataset);

    // we merge one partition each round
    for parts_left in tqdm!((0usize..9).rev()) {
        let mut generation = vec![];
        for basis in 1..=parts_left {
            for into in 0..basis {
                generation.push(best_mcm.merge(basis, into));
            }
        }
        let candidate = generation
            .iter()
            .max_by(|&a, &b| a.log_e(&dataset).total_cmp(&b.log_e(&dataset)))
            .unwrap();
        // if new >= old
        if candidate.log_e(&dataset).total_cmp(&best_log_e).is_ge() {
            best_mcm = candidate.clone();
            best_log_e = best_mcm.log_e(&dataset);
            println!("{best_mcm}");
            println!("With log E: {best_log_e}");
        } else {
            break;
        }
    }

    let mcm = best_mcm;
    println!("Best MCM with log E {}", best_log_e);
    println!("{}", mcm);
}
