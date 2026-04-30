use miette::Result;
use std::path::Path;

use mcm_finder_lib::{AnnealingStarter, AnnealingTemperature, SimulatedAnnealingSearcher, Solver};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = SimulatedAnnealingSearcher::from_file(filepath)?
        .set_temperature(AnnealingTemperature::new(10_000.0, 5.0, 0.0003).then(0.001, 0.00001))
        // .set_starter(AnnealingStarter::Trivial);
        .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
