use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{GreedySearcher, Solver};

fn main() -> Result<()> {
    // let dataset = Dataset::read_from_file(Path::new(
    //     "mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat",
    // ))
    // .unwrap();
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = GreedySearcher::from_file(filepath)?;
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
