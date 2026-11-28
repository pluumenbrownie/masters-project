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

pub(crate) fn parallel_tempering() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 10])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_shorter() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(8 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 10])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_-1_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_-1_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_-1_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_-1_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_less_short() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(9 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 10])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_-0_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_-0_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_-0_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_-0_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_longer() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(11 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 2])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+2_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+2_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+2_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+2_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_tabu() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_tabu";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 3])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_t_2048_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_t_2048_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_t_2048_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_tabu_list(Some(2048))
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_t_2048_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_hotter() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_hotter";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 3])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_t_10_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_t_10_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_t_10_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(10.0)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_t_10_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_less_pools() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_less_pools";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 3])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_t_6_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_t_6_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_t_6_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_pools(6)
                .set_arbitrary_swaps(false);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_t_6_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_random_swaps() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_random_swaps";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 3])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(true);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_t_1_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_hotter_peak() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_hotter_peak";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 5])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(3000.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(true);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_t_1_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_colder_peak() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_colder_peak";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|nr| [nr; 5])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles))| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_t_1_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(300.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_arbitrary_swaps(true);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_t_1_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}

pub(crate) fn parallel_tempering_boost_variable() -> Result<()> {
    const FOLDER: &str = "results/parallel_tempering_boost_variable";
    let dataset = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");

    let iterator: Vec<_> = (0..6)
        .into_iter()
        .map(|x| (2usize.pow(5 + x), 2usize.pow(10 - x)))
        .flat_map(|t| {
            [
                (ParallelTemperatureCurve::Linear, "l"),
                (ParallelTemperatureCurve::InverseLinear, "i"),
                (ParallelTemperatureCurve::Geometric, "g"),
                (ParallelTemperatureCurve::Exponential, "e"),
            ]
            .into_iter()
            .map(move |c| (c, t))
        })
        .flat_map(|s| [16usize, 128].into_iter().map(move |c| (s.0, s.1, c)))
        .flat_map(|nr| [nr; 5])
        .collect();

    par_tqdm!(iterator.into_par_iter()).try_for_each(
        |((curve, curve_name), (per_shuffle, shuffles), boost_variable)| -> Result<()> {
            let id = gen_id();

            let (collector, sender) = CollectorBuilder::<SolverEvent>::new()
                .with_handler(mcm_finder_lib::logger::ParTempHandler::with_paths(
                    format!("./{FOLDER}/log_e_{curve_name}_{per_shuffle}_+1_b_{boost_variable}_{id}.csv").into(),
                    format!("./{FOLDER}/swap_history_{curve_name}_{per_shuffle}_+1_b_{boost_variable}_{id}.csv").into(),
                    format!("./{FOLDER}/mcm_id_{curve_name}_{per_shuffle}_+1_b_{boost_variable}_{id}.csv").into(),
                ))
                .build();
            std::thread::spawn(|| collector.run());

            let solver = ParallelTemperingSolver::<VecDataset>::from_file(dataset)?
                .set_sender(sender)
                .set_starter(Trivial)
                .set_temperature_curve(curve)
                .set_max_temperature(900.0)
                .set_min_temperature(0.1)
                .set_steps_per_shuffle(per_shuffle)
                .set_count_non_cached(false)
                .set_shuffles(shuffles)
                .set_swap_acceptance_booster(boost_variable as f64);
            let result = solver.solve();
            result
                .save_json(format!("./{FOLDER}/result_{curve_name}_{per_shuffle}_+1_b_{boost_variable}_{id}.json").into())
                .into_diagnostic()?;
            Ok(())
        },
    )?;

    Ok(())
}
