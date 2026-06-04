use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{ExhaustiveSolver, Solver};

fn main() -> Result<()> {
    let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");

    let solver = ExhaustiveSolver::from_file(filepath);
    let result = solver?.solve();
    println!("{result}");
    println!();
    let matrix = result.mcm.to_matrix();
    for row in matrix {
        println!("{}", row);
    }
    Ok(())
}
