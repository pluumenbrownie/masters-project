use std::path::Path;

use miette::Result;

use mcm_finder_lib::{
    logger::get_collector,
    solvers::{
        ConstructiveStrategy, GreedySolver,
        InitialSolver::{Construct, ConstructBeam, Merge, MergeBeam, None},
        Refinements, Solver,
    },
};

fn main() -> Result<()> {
    let (collector, tx) = get_collector();
    std::thread::spawn(|| collector.run());

    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST14.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST22.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST28.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = GreedySolver::from_file(filepath)?
        .set_sender(tx)
        .set_initial_solver(Merge)
        // .set_refinement_sequence(vec![Refinements::LocalBeam { beam_size: 5 }])
        .set_refinement_sequence(vec![
            // Refinements::ChooseN {
            //     n: 4,
            //     max_fails: 500,
            // },
            // Refinements::Local,
            // Refinements::Exhaustive,
            Refinements::Tabu {
                steps: 1000,
                size: 100,
            },
            Refinements::Local,
        ]);
    // let solver = GreedySearcher::from_file(filepath)?.continue_after_minimum();
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
