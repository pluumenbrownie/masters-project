use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{
    AnnealingStarter, ParallelTemperatureCurve, ParallelTemperingSearcher, Solver,
};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = ParallelTemperingSearcher::from_file(filepath)?
        .set_starter(AnnealingStarter::Trivial)
        .set_temperature_curve(ParallelTemperatureCurve::Exponential);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
