use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{GreedySearcher, Solver};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = GreedySearcher::from_file(filepath)?.lookahead(0);
    // let solver = GreedySearcher::from_file(filepath)?.continue_after_minimum();
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
