use rand::seq::IndexedRandom;
use std::{any::type_name, hash::BuildHasher, path::Path, time::Duration};

use criterion::{
    BatchSize, BenchmarkId, Criterion, SamplingMode::Linear, criterion_group, criterion_main,
};

use fixedbitset::FixedBitSet;
use mcm_finder_lib::dataset::{
    AhashState, DataContainter, DataMap, DataVec, Dataset, DefaultState, EndsCachedDataset,
    EndsCachedVecDataset, FxState, RapidStateFast, RapidStateQuality, SimpleDataset, VecDataset,
};

fn criterion_benchmark(c: &mut Criterion) {
    {
        let mut c = c.benchmark_group("SimpleDataset");
        c = simple_dataset_benches::<DataVec, DefaultState>(c);
        c = simple_dataset_benches::<DataVec, AhashState>(c);
        c = simple_dataset_benches::<DataVec, RapidStateFast>(c);
        c = simple_dataset_benches::<DataVec, RapidStateQuality>(c);
        c = simple_dataset_benches::<DataVec, FxState>(c);
        c = simple_dataset_benches::<DataMap<DefaultState>, DefaultState>(c);
        c = simple_dataset_benches::<DataMap<AhashState>, AhashState>(c);
        c = simple_dataset_benches::<DataMap<RapidStateFast>, RapidStateFast>(c);
        c = simple_dataset_benches::<DataMap<RapidStateQuality>, RapidStateQuality>(c);
        c = simple_dataset_benches::<DataMap<FxState>, FxState>(c);
    }
    {
        let mut c = c.benchmark_group("Basic EndsCachedDataset");
        c = ends_cached_dataset_benches::<DataVec, DefaultState>(c);
        c = ends_cached_dataset_benches::<DataVec, AhashState>(c);
        c = ends_cached_dataset_benches::<DataVec, RapidStateFast>(c);
        c = ends_cached_dataset_benches::<DataVec, RapidStateQuality>(c);
        c = ends_cached_dataset_benches::<DataVec, FxState>(c);
        c = ends_cached_dataset_benches::<DataMap<DefaultState>, DefaultState>(c);
        c = ends_cached_dataset_benches::<DataMap<AhashState>, AhashState>(c);
        c = ends_cached_dataset_benches::<DataMap<RapidStateFast>, RapidStateFast>(c);
        c = ends_cached_dataset_benches::<DataMap<RapidStateQuality>, RapidStateQuality>(c);
        c = ends_cached_dataset_benches::<DataMap<FxState>, FxState>(c);
    }
    {
        let dataset = VecDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap();
        let mut group = c.benchmark_group("VecDataset Log E MNIST11 random inputs");
        for nr_of_bits in 0usize..122 {
            group.bench_with_input(
                BenchmarkId::from_parameter(nr_of_bits),
                &nr_of_bits,
                |b, &nr_of_bits| {
                    b.iter_batched(
                        || random_fixedbitset(122, nr_of_bits),
                        |icc| dataset.log_e(&icc),
                        BatchSize::SmallInput,
                    );
                },
            );
        }
        group.finish();
    }
    {
        let dataset =
            EndsCachedVecDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap();
        let mut group = c.benchmark_group("EndsCachedDataset Log E MNIST11 random inputs");
        for nr_of_bits in 0usize..122 {
            group.bench_with_input(
                BenchmarkId::from_parameter(nr_of_bits),
                &nr_of_bits,
                |b, &nr_of_bits| {
                    b.iter_batched(
                        || random_fixedbitset(122, nr_of_bits),
                        |icc| dataset.log_e(&icc),
                        BatchSize::SmallInput,
                    );
                },
            );
        }
        group.finish();
    }
}

fn ends_cached_dataset_benches<C: DataContainter<S>, S: BuildHasher + Default>(
    mut c: criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) -> criterion::BenchmarkGroup<'_, criterion::measurement::WallTime> {
    c.sample_size(10).measurement_time(Duration::from_secs(30));

    c.bench_function(format!("{}: Loading MNIST11", type_name::<C>()), |b| {
        b.iter(|| {
            EndsCachedDataset::<C, S>::read_from_file(Path::new("tests/data/MNIST11.sorted"))
                .unwrap()
        })
    });
    c.bench_function(format!("{}: Loading Big5", type_name::<C>()), |b| {
        b.iter(|| {
            EndsCachedDataset::<C, S>::read_from_file(Path::new("tests/data/Big5PT.sorted"))
                .unwrap()
        })
    });

    let dataset =
        SimpleDataset::<C, S>::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap();
    let icc = FixedBitSet::with_capacity_and_blocks(
        50,
        [0b0000000001000000000100000001000000000100000000010],
    );
    c.bench_function(format!("{}: Calculate Log E Big5", type_name::<C>()), |b| {
        b.iter(|| dataset.log_e(&icc))
    });
    c
}

fn simple_dataset_benches<C: DataContainter<S>, S: BuildHasher + Default>(
    mut c: criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
) -> criterion::BenchmarkGroup<'_, criterion::measurement::WallTime> {
    c.sample_size(10).measurement_time(Duration::from_secs(30));

    c.bench_function(
        format!(
            "{} + {}: Loading MNIST11",
            type_name::<C>(),
            type_name::<S>()
        ),
        |b| {
            b.iter(|| {
                SimpleDataset::<C, S>::read_from_file(Path::new("tests/data/MNIST11.sorted"))
                    .unwrap()
            })
        },
    );
    c.bench_function(
        format!("{} + {}: Loading Big5", type_name::<C>(), type_name::<S>()),
        |b| {
            b.iter(|| {
                SimpleDataset::<C, S>::read_from_file(Path::new("tests/data/Big5PT.sorted"))
                    .unwrap()
            })
        },
    );

    let dataset =
        SimpleDataset::<C, S>::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap();
    let icc = FixedBitSet::with_capacity_and_blocks(
        50,
        [0b0000000001000000000100000001000000000100000000010],
    );
    c.bench_function(format!("{}: Calculate Log E Big5", type_name::<C>()), |b| {
        b.iter(|| dataset.log_e(&icc))
    });
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
