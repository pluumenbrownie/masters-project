use miette::Result;
use std::path::Path;

use mcm_finder_lib::{GreedySearcher, Solver};

fn main() -> Result<()> {
    // let dataset = Dataset::read_from_file(Path::new(
    //     "mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat",
    // ))
    // .unwrap();
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let dataset =
    //     Dataset::read_from_file(Path::new("mcm-finder-lib/tests/data/Big5PT.sorted")).unwrap();

    let solver = GreedySearcher::from_file(filepath)?;
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
