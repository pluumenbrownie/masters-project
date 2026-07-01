use std::{
    any::{type_name, type_name_of_val},
    hash::BuildHasher,
    path::Path,
    time::Duration,
};

use criterion::{
    BatchSize, BenchmarkId, Criterion, SamplingMode::Linear, criterion_group, criterion_main,
};
use fixedbitset::FixedBitSet;
use rand::seq::IndexedRandom;

use mcm_finder_lib::dataset::{
    AhashState, Dataset, DefaultState, FxState, RapidStateFast, RapidStateQuality,
    datacontainer::{DataContainer, DataMap, DataVec},
    ends_cached::EndsCachedDataset,
    simple::SimpleDataset,
};

fn criterion_benchmark(c: &mut Criterion) {
    {
        // let datapaths = [DataPath::Scotus];
        let datapaths = [DataPath::Mnist11, DataPath::Big5];
        for data in datapaths {
            let mut c = c.benchmark_group(format!("SimpleDataset - Load {}", data.name()));
            c = simple_dataset_loading_bench::<DataVec<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataVec<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataVec<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataVec<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataVec<FxState>, FxState>(
                data,
                "FxState".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataMap<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataMap<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataMap<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataMap<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = simple_dataset_loading_bench::<DataMap<FxState>, FxState>(
                data,
                "FxState".into(),
                c,
            );
            c.finish();
        }
    }
    {
        // let datapaths = [DataPath::Scotus];
        let datapaths = [DataPath::Mnist11, DataPath::Big5];
        for data in datapaths {
            let mut c = c.benchmark_group(format!(
                "SimpleDataset - Calculate Log E of {}",
                data.name()
            ));
            c = simple_dataset_log_e_bench::<DataVec<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataVec<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataVec<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataVec<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataVec<FxState>, FxState>(data, "FxState".into(), c);
            c = simple_dataset_log_e_bench::<DataMap<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataMap<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataMap<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataMap<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = simple_dataset_log_e_bench::<DataMap<FxState>, FxState>(data, "FxState".into(), c);
            c.finish();
        }
    }
    {
        // let datapaths = [DataPath::Scotus];
        let datapaths = [DataPath::Mnist11, DataPath::Big5];
        for data in datapaths {
            let mut c = c.benchmark_group(format!("EndsCachedDataset - Load {}", data.name()));
            c = ends_cached_dataset_loading_bench::<DataVec<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataVec<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataVec<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataVec<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataVec<FxState>, FxState>(
                data,
                "FxState".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataMap<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataMap<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataMap<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataMap<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = ends_cached_dataset_loading_bench::<DataMap<FxState>, FxState>(
                data,
                "FxState".into(),
                c,
            );
            c.finish();
        }
    }
    {
        // let datapaths = [DataPath::Scotus];
        let datapaths = [DataPath::Mnist11, DataPath::Big5];
        for data in datapaths {
            let mut c = c.benchmark_group(format!(
                "EndsCachedDataset - Calculate Log E of {}",
                data.name()
            ));
            c = ends_cached_dataset_log_e_bench::<DataVec<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataVec<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataVec<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataVec<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataVec<FxState>, FxState>(
                data,
                "FxState".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataMap<DefaultState>, DefaultState>(
                data,
                "DefaultState".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataMap<AhashState>, AhashState>(
                data,
                "AhashState".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataMap<RapidStateFast>, RapidStateFast>(
                data,
                "RapidStateFast".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataMap<RapidStateQuality>, RapidStateQuality>(
                data,
                "RapidStateQuality".into(),
                c,
            );
            c = ends_cached_dataset_log_e_bench::<DataMap<FxState>, FxState>(
                data,
                "FxState".into(),
                c,
            );
            c.finish();
        }
    }
    // {
    //     let dataset = VecDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap();
    // let mut group = c.benchmark_group("VecDataset Log E MNIST11 random inputs");
    //     for nr_of_bits in 0usize..122 {
    //         group.bench_with_input(
    //             BenchmarkId::from_parameter(nr_of_bits),
    //             &nr_of_bits,
    //             |b, &nr_of_bits| {
    //                 b.iter_batched(
    //                     || random_fixedbitset(122, nr_of_bits),
    //                     |icc| dataset.log_e(&icc),
    //                     BatchSize::SmallInput,
    //                 );
    //             },
    //         );
    //     }
    //     group.finish();
    // }
    // {
    //     let dataset =
    //         EndsCachedVecDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap();
    //     let mut group = c.benchmark_group("EndsCachedDataset Log E MNIST11 random inputs");
    //     for nr_of_bits in 0usize..122 {
    //         group.bench_with_input(
    //             BenchmarkId::from_parameter(nr_of_bits),
    //             &nr_of_bits,
    //             |b, &nr_of_bits| {
    //                 b.iter_batched(
    //                     || random_fixedbitset(122, nr_of_bits),
    //                     |icc| dataset.log_e(&icc),
    //                     BatchSize::SmallInput,
    //                 );
    //             },
    //         );
    //     }
    //     group.finish();
    // }
}

#[derive(Clone, Copy)]
enum DataPath {
    Scotus,
    Mnist11,
    Big5,
}

impl DataPath {
    fn name(&self) -> String {
        match self {
            DataPath::Scotus => String::from("SCOTUS"),
            DataPath::Mnist11 => String::from("MNIST11"),
            DataPath::Big5 => String::from("Big5"),
        }
    }

    fn get_filepath(&self) -> &Path {
        match self {
            DataPath::Scotus => Path::new("tests/data/SCOTUS_n9_N895_Data.dat"),
            DataPath::Mnist11 => Path::new("tests/data/MNIST11.sorted"),
            DataPath::Big5 => Path::new("tests/data/Big5PT.sorted"),
        }
    }
}

fn simple_dataset_loading_bench<'a, C: DataContainer<S>, S: BuildHasher + Default>(
    datapath: DataPath,
    hasher_name: String,
    mut c: criterion::BenchmarkGroup<'a, criterion::measurement::WallTime>,
) -> criterion::BenchmarkGroup<'a, criterion::measurement::WallTime> {
    c.sample_size(10).measurement_time(Duration::from_secs(30));

    c.bench_with_input(
        BenchmarkId::new(format!("Loading {}", C::name()), hasher_name),
        &(),
        |b, _| b.iter(|| SimpleDataset::<C, S>::read_from_file(datapath.get_filepath()).unwrap()),
    );

    c
}

fn simple_dataset_log_e_bench<C: DataContainer<S>, S: BuildHasher + Default>(
    datapath: DataPath,
    hasher_name: String,
    mut c: criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) -> criterion::BenchmarkGroup<'_, criterion::measurement::WallTime> {
    c.sample_size(10).measurement_time(Duration::from_secs(30));

    let dataset = SimpleDataset::<C, S>::read_from_file(datapath.get_filepath()).unwrap();
    let icc = match datapath {
        DataPath::Scotus => FixedBitSet::with_capacity_and_blocks(9, [0b010001001]),
        _ => FixedBitSet::with_capacity_and_blocks(
            50,
            [0b0000000001000000000100000001000000000100000000010],
        ),
    };

    c.bench_with_input(
        BenchmarkId::new(format!("Calculate Log E with {}", C::name()), hasher_name),
        &(),
        |b, _| b.iter(|| dataset.log_e(&icc)),
    );
    c
}

fn ends_cached_dataset_loading_bench<C: DataContainer<S>, S: BuildHasher + Default>(
    datapath: DataPath,
    hasher_name: String,
    mut c: criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) -> criterion::BenchmarkGroup<'_, criterion::measurement::WallTime> {
    c.sample_size(10).measurement_time(Duration::from_secs(30));

    c.bench_with_input(
        BenchmarkId::new(format!("Loading {}", C::name()), hasher_name),
        &(),
        |b, _| {
            b.iter(|| EndsCachedDataset::<C, S>::read_from_file(datapath.get_filepath()).unwrap())
        },
    );
    c
}

fn ends_cached_dataset_log_e_bench<C: DataContainer<S>, S: BuildHasher + Default>(
    datapath: DataPath,
    hasher_name: String,
    mut c: criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) -> criterion::BenchmarkGroup<'_, criterion::measurement::WallTime> {
    c.sample_size(10).measurement_time(Duration::from_secs(30));

    let dataset = EndsCachedDataset::<C, S>::read_from_file(datapath.get_filepath()).unwrap();
    let icc = FixedBitSet::with_capacity_and_blocks(
        50,
        [0b0000000001000000000100000001000000000100000000010],
    );
    c.bench_with_input(
        BenchmarkId::new(format!("Calculate Log E with {}", C::name()), hasher_name),
        &(),
        |b, _| b.iter(|| dataset.log_e(&icc)),
    );
    c
}

fn random_fixedbitset(length: usize, nr_of_bits: usize) -> FixedBitSet {
    let mut icc = FixedBitSet::with_capacity(length);
    icc.set_range(.., false);
    let mut rng = rand::rng();
    let numbers: Vec<_> = (0..length).collect();
    for nr in numbers.sample(&mut rng, nr_of_bits) {
        icc.set(*nr, true);
    }
    icc
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
