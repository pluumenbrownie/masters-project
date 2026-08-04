use std::path::Path;

use miette::Result;

use mcm_finder_lib::{
    dataset::MapDataset,
    solvers::{
        AntColonyOptimization, ConstructiveStrategy, GreedySolver,
        InitialSolver::{Construct, ConstructBeam, Merge, MergeBeam, None},
        Refinements, Solver,
    },
};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST14.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST22.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST28.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver: AntColonyOptimization<MapDataset> = AntColonyOptimization::from_file(filepath)?;
    // let solver = GreedySearcher::from_file(filepath)?.continue_after_minimum();
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
