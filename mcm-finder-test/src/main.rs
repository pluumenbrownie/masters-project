use std::path::Path;

use annolog::CollectorBuilder;
use mcm_finder_lib::{
    dataset::VecDataset,
    logger::SolverEvent,
    solvers::{
        AnnealingStarter::Trivial, AnnealingTemperature, GreedySolver, HillClimberSolver,
        InitialSolver::Merge, Refinements, SimulatedAnnealingSolver, Solver,
    },
};

use clap::{Parser, Subcommand};
use kdam::par_tqdm;
use miette::{IntoDiagnostic, Result};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Greedy,
    Annealing,
    AnnealingTiming,
    Lahc,
}

fn main() -> Result<()> {
    let args = Cli::try_parse().unwrap_or_else(|error| error.exit());

    match args.command {
        Commands::Greedy => greedy_test()?,
        Commands::Annealing => annealing_length()?,
        Commands::AnnealingTiming => annealing_timing()?,
        Commands::Lahc => lahc_hist_size()?,
    };

    Ok(())
}

fn annealing_length() -> Result<()> {
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

fn lahc_hist_size() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = [
        // 50, 75, 100, 125, 150, 175,
        // 20, 30, 40,
        1, 2, 5, 10,
    ]
    .into_iter()
    .flat_map(|nr| [nr; 25])
    .enumerate()
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(nr, hist_size)| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::HillClimbingHandler::with_paths(
                format!("./results/lahc_hist_size/log_e_{nr}_{hist_size}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = HillClimberSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_starter(Trivial)
            .set_max_steps(500_000)
            .set_stagnation_steps(100_000)
            .set_history_size(hist_size);

        let result = solver.solve();
        result
            .save_json(format!("./results/lahc_hist_size/result_{nr}_{hist_size}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

fn annealing_timing() -> Result<()> {
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

fn greedy_test() -> Result<()> {
    let datasets = [
        Path::new("mcm-finder-lib/tests/data/Immobilized184neur/worm3.dat"),
        Path::new("mcm-finder-lib/tests/data/MNIST11.sorted"),
    ];

    for (nr, dataset) in datasets.iter().enumerate() {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                format!("log_e_greedy_{nr}").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = GreedySolver::from_file(dataset)?
            .set_sender(sender)
            .set_initial_solver(Merge)
            .set_refinement_sequence(vec![
                Refinements::Tabu {
                    steps: 1000,
                    size: 100,
                },
                Refinements::Local,
            ]);
    }

    Ok(())
}
