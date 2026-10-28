use std::path::Path;

use annolog::CollectorBuilder;
use mcm_finder_lib::{
    dataset::VecDataset,
    logger::SolverEvent,
    solvers::{
        AnnealingStarter::Trivial, ParallelTemperatureCurve, ParallelTemperingSolver, Solver,
    },
};

use kdam::{par_tqdm, tqdm};
use miette::{IntoDiagnostic, Result};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::gen_id;

pub(crate) fn parallel_timing() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_timing";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    for nr in tqdm!(0..100) {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                format!("./{FOLDER}/log_e_{nr}.csv").into(),
                format!("./{FOLDER}/swap_history_{nr}.csv").into(),
                format!("./{FOLDER}/mcm_id_{nr}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_starter(Trivial)
            .set_temperature_curve(ParallelTemperatureCurve::Linear)
            .set_max_temperature(900.0)
            .set_min_temperature(0.1)
            .set_steps_per_shuffle(50)
            .set_count_non_cached(false)
            .set_shuffles(1000)
            .set_arbitrary_swaps(false);
        let result = solver.solve();
        result
            .save_json(format!("./{FOLDER}/result_{nr}.json").into())
            .into_diagnostic()?;
    }

    Ok(())
}
