use std::path::Path;
use miette::Result;
use mcm_finder_lib::{
    solvers::{
        ConstructiveStrategy, GreedySolver,
        InitialSolver::{Construct, ConstructBeam, Merge, MergeBeam, None},
        Refinements, Solver,
    },
};

fn main() -> Result<()> {
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    let solver = GreedySolver::from_file(filepath)?
        .set_initial_solver(Merge)
        .set_refinement_sequence(vec![
            Refinements::Tabu {
                steps: 1000,
                size: 100,
            },
            Refinements::Local,
        ]);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}