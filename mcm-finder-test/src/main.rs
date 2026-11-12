use std::path::Path;

use annolog::CollectorBuilder;
use mcm_finder_lib::{
    dataset::VecDataset,
    logger::SolverEvent,
    solvers::{
        AnnealingStarter::Trivial, GreedySolver, HillClimberSolver, InitialSolver, Refinements,
        Solver,
    },
};

use clap::{Parser, Subcommand};
use kdam::{par_tqdm, tqdm};
use miette::{IntoDiagnostic, Result};

use rand::{RngExt, distr::Alphanumeric};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

mod annealing;
mod parallel;

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Greedy,
    TabuDivergence,
    TabuDivergenceLower,
    PureTabu,
    MergeBeam,
    ConstructBeam,
    Annealing,
    AnnealingTiming,
    AnnealingShapes,
    Adaptive,
    ParallelTiming,
    ParallelTempering,
    ParallelTemperingShorter,
    ParallelTemperingLessShort,
    ParallelTemperingLonger,
    Lahc,
}

fn main() -> Result<()> {
    let args = Cli::try_parse().unwrap_or_else(|error| error.exit());

    match args.command {
        Commands::Greedy => greedy_test()?,
        Commands::TabuDivergence => tabu_test()?,
        Commands::TabuDivergenceLower => tabu_test_2()?,
        Commands::PureTabu => pure_tabu()?,
        Commands::MergeBeam => merge_beam()?,
        Commands::ConstructBeam => construct_beam_timing()?,
        Commands::Annealing => annealing::annealing_length()?,
        Commands::AnnealingTiming => annealing::annealing_timing()?,
        Commands::AnnealingShapes => annealing::annealing_shapes()?,
        Commands::Adaptive => annealing::adaptive_annealing()?,
        Commands::ParallelTiming => parallel::parallel_timing()?,
        Commands::ParallelTempering => parallel::parallel_tempering()?,
        Commands::ParallelTemperingShorter => parallel::parallel_tempering_shorter()?,
        Commands::ParallelTemperingLessShort => parallel::parallel_tempering_less_short()?,
        Commands::ParallelTemperingLonger => parallel::parallel_tempering_longer()?,
        Commands::Lahc => lahc_hist_size()?,
    };

    Ok(())
}

pub fn gen_id() -> String {
    let id: String = rand::rng()
        .sample_iter(Alphanumeric)
        .take(8)
        .map(char::from)
        .collect();
    id
}

fn lahc_hist_size() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = [
        // 1, 2, 10,
        // 20,
        // 50,
        100,
        // 50, 200
    ]
    .into_iter()
    .flat_map(|nr| [nr; 100])
    .enumerate()
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(nr, hist_size)| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::HillClimbingHandler::with_paths(
                format!("./results/lahc_hist_size/log_e_a_{nr}_{hist_size}.csv").into(),
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
            .save_json(format!("./results/lahc_hist_size/result_a_{nr}_{hist_size}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

fn greedy_test() -> Result<()> {
    let datasets = [
        // (
        //     Path::new("mcm-finder-lib/tests/data/Immobilized184neur/worm3.dat"),
        //     "worm3",
        // ),
        (
            Path::new("mcm-finder-lib/tests/data/MNIST11.sorted"),
            "mnist11",
        ),
    ];

    let settings = [
        (InitialSolver::Merge, "m_1"),
        (InitialSolver::MergeBeam { beam_size: 2 }, "m_2"),
        (InitialSolver::MergeBeam { beam_size: 5 }, "m_5"),
        (InitialSolver::MergeBeam { beam_size: 20 }, "m_20"),
        (InitialSolver::MergeBeam { beam_size: 100 }, "m_100"),
        (InitialSolver::Construct, "c_1"),
        (InitialSolver::ConstructBeam { beam_size: 2 }, "c_2"),
        (InitialSolver::ConstructBeam { beam_size: 5 }, "c_5"),
        (InitialSolver::ConstructBeam { beam_size: 20 }, "c_20"),
        (InitialSolver::ConstructBeam { beam_size: 100 }, "c_100"),
        (InitialSolver::None, "n"),
    ];

    let refinements = [
        (
            vec![
                Refinements::Tabu {
                    steps: 10_000,
                    size: 100,
                },
                Refinements::Local,
            ],
            "t_100",
        ),
        (
            vec![
                Refinements::Tabu {
                    steps: 10_000,
                    size: 50,
                },
                Refinements::Local,
            ],
            "t_50",
        ),
        (
            vec![
                Refinements::Tabu {
                    steps: 10_000,
                    size: 25,
                },
                Refinements::Local,
            ],
            "t_25",
        ),
        (
            vec![
                Refinements::Tabu {
                    steps: 10_000,
                    size: 10,
                },
                Refinements::Local,
            ],
            "t_10",
        ),
        (vec![], "nr"),
        (vec![Refinements::Local], "l_1"),
        (vec![Refinements::LocalBeam { beam_size: 2 }], "l_2"),
        (vec![Refinements::LocalBeam { beam_size: 5 }], "l_5"),
        (vec![Refinements::LocalBeam { beam_size: 20 }], "l_20"),
        (vec![Refinements::LocalBeam { beam_size: 100 }], "l_100"),
    ];

    let iterator: Vec<_> = settings
        .into_iter()
        .flat_map(|s| datasets.iter().map(move |d| (d, s.clone())))
        .flat_map(|(d, s)| refinements.iter().map(move |r| (d, s.clone(), r)))
        .filter(|(_, s, r)| !(s.0 == InitialSolver::None && r.0.is_empty()))
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((dataset, dataset_name), (strategy, name), (refinements, ref_name))| -> Result<()> {
            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                    format!("./results/greedy/log_e_{dataset_name}_{name}_{ref_name}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = GreedySolver::from_file(dataset)?
                .set_sender(sender)
                .set_silent(true)
                .set_initial_solver(strategy)
                .set_refinement_sequence(refinements.clone());
            let result = solver.solve();
            result
                .save_json(
                    format!("./results/greedy/result_{dataset_name}_{name}_{ref_name}.json").into(),
                )
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

fn tabu_test() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    const FOLDER: &str = "results/tabu_divergence";

    let iterator: Vec<_> = [
        90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
        110,
    ]
    .into_iter()
    .map(|q| {
        (
            vec![
                Refinements::Tabu {
                    steps: 10_000,
                    size: q,
                },
                Refinements::Local,
            ],
            format!("t_{q}"),
        )
    })
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(refinements, name)| -> Result<()> {
        // for (nr, dataset) in datasets.iter().enumerate() {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                format!("./{FOLDER}/log_e_{name}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = GreedySolver::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_initial_solver(InitialSolver::ConstructBeam { beam_size: 2 })
            .set_refinement_sequence(refinements.clone());
        let result = solver.solve();
        result
            .save_json(format!("./{FOLDER}/result_{name}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

fn tabu_test_2() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    const FOLDER: &str = "results/tabu_divergence_lower";

    let iterator: Vec<_> = [
        // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    ]
    .into_iter()
    .map(|q| {
        (
            vec![
                Refinements::Tabu {
                    steps: 10_000,
                    size: q,
                },
                Refinements::Local,
            ],
            format!("t_{q}"),
        )
    })
    .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(refinements, name)| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                format!("./{FOLDER}/log_e_{name}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = GreedySolver::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_initial_solver(InitialSolver::ConstructBeam { beam_size: 2 })
            .set_refinement_sequence(refinements.clone());
        let result = solver.solve();
        result
            .save_json(format!("./{FOLDER}/result_{name}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

fn pure_tabu() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    const FOLDER: &str = "results/pure_tabu";

    let iterator: Vec<_> = (1..=100)
        .into_iter()
        .map(|q| {
            (
                vec![
                    Refinements::Tabu {
                        steps: 10_000,
                        size: q,
                    },
                    Refinements::Local,
                ],
                format!("t_{q}"),
            )
        })
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(refinements, name)| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                format!("./{FOLDER}/log_e_{name}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = GreedySolver::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_initial_solver(InitialSolver::None)
            .set_refinement_sequence(refinements.clone());
        let result = solver.solve();
        result
            .save_json(format!("./{FOLDER}/result_{name}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

fn merge_beam() -> Result<()> {
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    const FOLDER: &str = "results/merge_beam";

    let iterator: Vec<_> = (1..=100)
        .into_iter()
        .map(|q| (InitialSolver::MergeBeam { beam_size: q }, format!("m_{q}")))
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(|(initial, name)| -> Result<()> {
        let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
            .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                format!("./{FOLDER}/log_e_{name}.csv").into(),
            ))
            .build();
        std::thread::spawn(|| collector.run());

        let solver = GreedySolver::from_file(dataset)?
            .set_sender(sender)
            .set_silent(true)
            .set_initial_solver(initial);
        let result = solver.solve();
        result
            .save_json(format!("./{FOLDER}/result_{name}.json").into())
            .into_diagnostic()?;
        Ok(())
    })?;

    Ok(())
}

fn construct_beam_timing() -> Result<()> {
    let datasets = [
        // (
        //     Path::new("mcm-finder-lib/tests/data/MNIST11.sorted"),
        //     "mnist11",
        // ),
        (
            Path::new("mcm-finder-lib/tests/data/MNIST14.sorted"),
            "mnist14",
        ),
        (
            Path::new("mcm-finder-lib/tests/data/MNIST22.sorted"),
            "mnist22",
        ),
    ];
    const FOLDER: &str = "results/construct_beam";

    let iterator: Vec<_> = (1..=100)
        .into_iter()
        .map(|q| {
            (
                InitialSolver::ConstructBeam { beam_size: q },
                format!("c_{q}"),
            )
        })
        .flat_map(|s| datasets.iter().map(move |d| (s.clone(), d)))
        .collect();

    tqdm!(iterator.into_iter()).try_for_each(
        |((initial, name), (dataset, data_name))| -> Result<()> {
            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::GreedyHandler::with_paths(
                    format!("./{FOLDER}/log_e_{name}_{data_name}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = GreedySolver::from_file(dataset)?
                .set_sender(sender)
                .set_silent(true)
                .set_initial_solver(initial);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{name}_{data_name}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}
