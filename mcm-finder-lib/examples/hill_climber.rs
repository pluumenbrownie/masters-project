use std::path::Path;

use miette::Result;

use mcm_finder_lib::{
    dataset::{
        AhashState, SimpleDataset,
        datacontainer::{DataMap, DataVec, DataVecSoA},
        simple::VecDataset,
    },
    solvers::{AnnealingStarter, HillClimberSolver, Solver},
};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver =
        HillClimberSolver::<SimpleDataset<DataVec<AhashState>, AhashState>>::from_file(filepath)?
            .set_starter(AnnealingStarter::Trivial)
            .set_max_steps(10_000_000)
            .set_stagnation_steps(100_000)
            .set_history_size(100);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
