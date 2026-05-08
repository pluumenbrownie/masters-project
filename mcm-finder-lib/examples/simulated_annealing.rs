use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{
    AnnealingStarter, SimulatedAnnealingSearcher, Solver, anneal_temps::AnnealingTemperature,
};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    // let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = SimulatedAnnealingSearcher::from_file(filepath)?
        .set_temperature(
            AnnealingTemperature::logarithmic(1_000_000.0, 1_000.0)
                // .then_constant(10_000)
                .then_exponential(0.001, 0.003),
            // AnnealingTemperature::logarithmic(1_000_000.0, 1.0),
            // .then_exponential(5.0, 0.0003)
            // .then_constant(10_000)
            // .then_exponential(0.001, 0.00001),
        )
        .set_starter(AnnealingStarter::Trivial);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
