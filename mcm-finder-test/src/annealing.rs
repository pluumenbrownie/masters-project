use std::path::Path;

use annolog::CollectorBuilder;
use mcm_finder_lib::{
    dataset::VecDataset,
    logger::SolverEvent,
    solvers::{
        AdaptiveAnnealingSolver, AnnealingStarter::Trivial, AnnealingTemperature,
        SimulatedAnnealingSolver, Solver,
    },
};

use kdam::par_tqdm;
use miette::{IntoDiagnostic, Result};
use rand::{RngExt, distr::Alphanumeric};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::gen_id;

pub(crate) fn annealing_length() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = [
        // 0.002, 0.0015, 0.001, 0.0007, 0.0005, 0.0003,
        0.0002, 0.00015, 0.0001,
    ]
    .into_iter()
    .flat_map(|nr| [nr; 25])
    .enumerate()
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(nr, per_step)| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::AnnealingHandler::with_paths(
                format!("./results/simulated_annealing_longer/log_e_{nr}_{per_step:.4}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = SimulatedAnnealingSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_starter(Trivial)
            .set_temperature(
                AnnealingTemperature::logarithmic(1_000_000.0, 1_000.0)
                    .then_constant(10_000)
                    .then_exponential(0.0001, per_step),
            );
        let result = solver.solve();
        result
            .save_json(
                format!("./results/simulated_annealing_longer/result_{nr}_{per_step:.4}.json")
                    .into(),
            )
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

pub(crate) fn annealing_shapes() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = [
        (
            AnnealingTemperature::exponential(100_000.0, 0.0001, 0.0001),
            "exp_01",
        ),
        (
            AnnealingTemperature::exponential(100_000.0, 0.0001, 0.0002),
            "exp_02",
        ),
        (
            AnnealingTemperature::exponential(100_000.0, 0.0001, 0.0005),
            "exp_05",
        ),
        (
            AnnealingTemperature::exponential(100_000.0, 0.0001, 0.001),
            "exp_1",
        ),
        (
            AnnealingTemperature::exponential(100_000.0, 0.0001, 0.002),
            "exp_2",
        ),
        (
            AnnealingTemperature::exponential(100_000.0, 0.0001, 0.01),
            "exp_10",
        ),
        (
            AnnealingTemperature::logarithmic(1_000_000.0, 10.0),
            "log_10",
        ),
        (
            AnnealingTemperature::logarithmic(1_000_000.0, 5.0),
            "log_05",
        ),
        (
            AnnealingTemperature::logarithmic(1_000_000.0, 3.0),
            "log_03",
        ),
        (
            AnnealingTemperature::logarithmic(100_000.0, 10.0),
            "log_10_s1",
        ),
        (
            AnnealingTemperature::logarithmic(200_000.0, 10.0),
            "log_10_s2",
        ),
        (
            AnnealingTemperature::logarithmic(500_000.0, 10.0),
            "log_10_s5",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 1.0),
            "lin_1",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 2.0),
            "lin_2",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 5.0),
            "lin_5",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.01, 10.0),
            "lin_10_h100",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 1.0, 10.0),
            "lin_10_h10000",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 100.0, 10.0),
            "lin_10_h1000000",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 10.0),
            "lin_10",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 20.0),
            "lin_20",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 50.0),
            "lin_50",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 100.0),
            "lin_100",
        ),
        (
            AnnealingTemperature::linear(100_000.0, 0.0001, 200.0),
            "lin_200",
        ),
        (
            AnnealingTemperature::logarithmic(1_000_000.0, 1_000.0)
                .then_constant(10_000)
                .then_exponential(0.0001, 0.0002),
            "cus_02",
        ),
    ]
    .into_iter()
    .flat_map(|t| [t].into_iter().cycle().take(25))
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(schedule, name)| -> Result<()> {
        let id = gen_id();

        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::AnnealingHandler::with_paths(
                format!("./results/simulated_annealing_shapes/log_e_{name}_{id}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = SimulatedAnnealingSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_starter(Trivial)
            .set_temperature(schedule);
        let result = solver.solve();
        result
            .save_json(
                format!("./results/simulated_annealing_shapes/result_{name}_{id}.json").into(),
            )
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

pub(crate) fn annealing_timing() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    for nr in 0..100 {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::AnnealingHandler::with_paths(
                format!("./results/simulated_annealing/log_e_{nr}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = SimulatedAnnealingSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_starter(Trivial)
            .set_temperature(
                AnnealingTemperature::logarithmic(1_000_000.0, 1_000.0)
                    .then_constant(10_000)
                    .then_exponential(0.0001, 0.002),
            );
        let result = solver.solve();
        result
            .save_json(format!("./results/simulated_annealing/result_{nr}.json").into())
            .into_diagnostic()?;
    }

    Ok(())
}

pub(crate) fn adaptive_annealing() -> Result<()> {
    let folder = "results/adaptive_annealing";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = [
        // 0.002, 0.0015, 0.001, 0.0007, 0.0005,
        0.003, 0.004,
    ]
    .into_iter()
    .flat_map(|d| {
        [0.95, 0.90, 0.85, 0.8, 0.75]
            .into_iter()
            .map(move |r| (d, r))
    })
    .flat_map(|t| [t].into_iter().cycle().take(25))
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(decay, reset)| -> Result<()> {
        let id = gen_id();

        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::AnnealingHandler::with_paths(
                format!("./{folder}/log_e_{decay}_{reset}_{id}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = AdaptiveAnnealingSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_starter(Trivial)
            .set_delta_per_step(decay)
            .set_reset_fraction(reset);
        let result = solver.solve();
        result
            .save_json(format!("./{folder}/result_{decay}_{reset}_{id}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}
