use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{AnnealingStarter, ParallelTemperingSearcher, Solver};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = ParallelTemperingSearcher::from_file(filepath)?
        .set_starter(AnnealingStarter::Trivial)
        .set_temperatures(vec![
            2000.0, 1000.0, 500.0, 250.0, 150.0, 50.0, 25.0, 12.0, 8.0, 5.0, 3.0, 1.0, 0.1,
        ]);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
