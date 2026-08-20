use std::path::Path;

use miette::Result;

use mcm_finder_lib::{
    logger::get_collector,
    solvers::{AnnealingStarter, ParallelTemperatureCurve, ParallelTemperingSolver, Solver},
};

fn main() -> Result<()> {
    let (collector, tx) = get_collector();
    std::thread::spawn(|| collector.run());

    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = ParallelTemperingSolver::from_file(filepath)?
        .set_sender(tx)
        .set_starter(AnnealingStarter::Trivial)
        .set_temperature_curve(ParallelTemperatureCurve::Geometric)
        .set_max_temperature(500.0)
        .set_min_temperature(0.1)
        .set_steps_per_shuffle(50)
        .set_count_non_cached(false)
        .set_shuffles(1500)
        .set_arbitrary_swaps(false);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
