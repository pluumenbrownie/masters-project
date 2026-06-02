use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{
    ConstructiveSolver,
    ConstructiveStrategy::{
        BackToFront, FrontToBack, GreedyOrder, RandomOrder, VarianceFirst, VarianceLast,
    },
    Solver,
};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST14.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST22.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST28.sorted");
    let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = ConstructiveSolver::from_file(filepath)?.set_strategy(FrontToBack);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
