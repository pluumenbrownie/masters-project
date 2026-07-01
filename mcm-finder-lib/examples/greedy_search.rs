use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{
    GreedySolver,
    InitialSolver::{Construct, Merge, None},
    Refinements, Solver,
};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST14.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST22.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST28.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = GreedySolver::from_file(filepath)?
        .lookahead(0)
        .set_initial_solver(Construct)
        .set_refinement_sequence(vec![
            Refinements::ChooseN {
                n: 4,
                max_fails: 500,
            },
            Refinements::Local,
        ]);
    // let solver = GreedySearcher::from_file(filepath)?.continue_after_minimum();
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
