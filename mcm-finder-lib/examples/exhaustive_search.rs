use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{ExhaustiveSearcher, Solver};

fn main() -> Result<()> {
    let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");

    let solver = ExhaustiveSearcher::from_file(filepath);
    let result = solver?.solve();
    println!("{result}");
    Ok(())
}
