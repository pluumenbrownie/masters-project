use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use annolog::CollectorBuilder;
use mcm_finder_lib::{
    dataset::VecDataset,
    logger::SolverEvent,
    solvers::{
        AnnealingStarter::Trivial, AnnealingTemperature, GreedySolver, InitialSolver::Merge,
        Refinements, SimulatedAnnealingSolver, Solver,
    },
};

use clap::{Parser, Subcommand};
use miette::Result;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use ron::ser::PrettyConfig;

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
}

fn main() -> Result<()> {
    let args = Cli::try_parse().unwrap_or_else(|error| error.exit());

    match args.command {
        Commands::Greedy => greedy_test()?,
        Commands::Annealing => annealing_test()?,
        Commands::AnnealingTiming => annealing_timing()?,
    };

    Ok(())
}

fn annealing_test() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    (0..100).into_par_iter().try_for_each(|nr| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::AnnealingHandler::with_paths(
                format!("./results/simulated_annealing_longer/log_e_{nr}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = SimulatedAnnealingSolver::<VecDataset>::from_file(dataset)?
            .set_sender(sender)
            .set_starter(Trivial)
            .set_temperature(
                AnnealingTemperature::logarithmic(1_000_000.0, 1_000.0)
                    .then_constant(10_000)
                    .then_exponential(0.0001, 0.001),
            );
        let result = solver.solve();
        let mut file = File::create(format!(
            "./results/simulated_annealing_longer/result_{nr}.ron"
        ))
        .unwrap();
        file.write_all(
            &ron::ser::to_string_pretty(&result, PrettyConfig::default())
                .unwrap()
                .into_bytes(),
        )
        .unwrap();
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
            .set_starter(Trivial)
            .set_temperature(
                AnnealingTemperature::logarithmic(1_000_000.0, 1_000.0)
                    .then_constant(10_000)
                    .then_exponential(0.0001, 0.002),
            );
        let result = solver.solve();
        let mut file =
            File::create(format!("./results/simulated_annealing/result_{nr}.ron")).unwrap();
        file.write_all(
            &ron::ser::to_string_pretty(&result, PrettyConfig::default())
                .unwrap()
                .into_bytes(),
        )
        .unwrap();
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
