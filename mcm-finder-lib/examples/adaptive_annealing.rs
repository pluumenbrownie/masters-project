use std::path::Path;

use miette::Result;

use mcm_finder_lib::{
    logger::get_collector,
    solvers::{AdaptiveAnnealingSolver, AnnealingStarter, Solver},
};

fn main() -> Result<()> {
    let (collector, tx) = get_collector();
    std::thread::spawn(|| collector.run());

    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = AdaptiveAnnealingSolver::from_file(filepath)?
        .set_sender(tx)
        .set_starter(AnnealingStarter::Trivial);
    // .set_starter(AnnealingStarter::Single);
    let result = solver.solve();
    println!("{}", result);
    Ok(())
}
