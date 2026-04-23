use std::path::Path;

use mcm_finder_lib::{ExhaustiveSearcher, Solver};
use miette::Result;

fn main() -> Result<()> {
    let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");

    let solver = ExhaustiveSearcher::from_file(filepath);
    let result = solver?.solve();
    println!("{result}");
    Ok(())
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
