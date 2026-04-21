use std::path::Path;

use fixedbitset::FixedBitSet;
use kdam::tqdm;
use mcm_finder_lib::*;

fn main() {
    let dataset = dataset::Dataset::read_from_file(Path::new(
        "mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat",
    ))
    .unwrap();

    let mcms = generate_all_mcms(9);
    let mut best_mcm: Option<MinimallyComplexModel> = None;
    let mut best_log_e = f64::NEG_INFINITY;

    for mcm in tqdm!(mcms.into_iter()) {
        let log_e = mcm.log_e(&dataset);
        if log_e > best_log_e {
            best_log_e = log_e;
            best_mcm = Some(mcm);
        }
    }

    let mcm = best_mcm.unwrap();
    println!("Best MCM with log E {}", best_log_e);
    println!("{}", mcm);
}

fn generate_all_mcms(number: usize) -> Vec<MinimallyComplexModel> {
    let mut output: Vec<Vec<usize>> = vec![];
    let mut current = vec![0usize; number];

    add_to(&mut current, 0, &mut output);

    output
        .into_iter()
        .map(|v| {
            let mut partitions = vec![
                FixedBitSet::with_capacity_and_blocks(9, [0b000000000]);
                *v.iter().max().unwrap() + 1
            ];
            for (nr, &i) in v.iter().enumerate() {
                partitions[i].set(nr, true);
            }
            MinimallyComplexModel::new(partitions)
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

// fn permutations(number: usize) -> Vec<Vec<usize>> {
//     let mut output = HashSet::new();

//     let mut old_generation = vec![vec![1usize; number]];
//     output.insert(old_generation[0].clone());

//     // run number-1 times
//     for _ in 1..number {
//         let new_generation: Vec<Vec<usize>> = old_generation
//             .iter()
//             .flat_map(|g| {
//                 // find the indexes we want to
//                 let mut out = vec![];
//                 let mut g = g.clone();

//                 let mut state = SplitState::SearchingOne;

//                 for nr in (0..g.len()).rev() {
//                     match state {
//                         SplitState::SearchingOne => {
//                             if g[nr] == 1 {
//                                 g[nr] = 0;
//                                 state = SplitState::Adding;
//                             }
//                         }
//                         SplitState::Adding => {
//                             if nr == 0 || g[nr] < g[nr - 1] {
//                                 let mut new_g = g.clone();
//                                 new_g[nr] += 1;
//                                 out.push(new_g);
//                             }
//                         }
//                     };
//                 }
//                 out
//             })
//             .collect();
//         old_generation.clear();
//         for new in new_generation {
//             dbg!(output.insert(new.clone()));
//             old_generation.push(new);
//         }
//     }

//     let mut out_vec: Vec<Vec<usize>> = Vec::from_iter(output.into_iter());
//     out_vec.sort();
//     out_vec
// }

// enum SplitState {
//     SearchingOne,
//     Adding,
// }
