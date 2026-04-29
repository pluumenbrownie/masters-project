use miette::Result;
use std::path::Path;

use mcm_finder_lib::{
    AnnealingStarter, SimulatedAnnealingSearcher, Solver, SolverReport, mcm::MinimallyComplexModel,
};

fn main() -> Result<()> {
    // let dataset = Dataset::read_from_file(Path::new(
    //     "mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat",
    // ))
    // .unwrap();
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = SimulatedAnnealingSearcher::from_file(filepath)?
        .set_temperature(50_000_000.0, 10_000.0, 0.001)
        .set_starter(AnnealingStarter::Trivial);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
