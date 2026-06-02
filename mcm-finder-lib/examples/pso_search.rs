use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{PsoSolver, PsoStrategy::Rounding, Solver};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST14.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST22.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST28.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = PsoSolver::from_file(filepath)?
        .set_strategy(Rounding)
        .set_momentum_weight(0.79);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
