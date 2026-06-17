use rand::seq::IndexedRandom;
use std::{path::Path, time::Duration};

use criterion::{
    BatchSize, BenchmarkId, Criterion, SamplingMode::Linear, criterion_group, criterion_main,
};

use fixedbitset::FixedBitSet;
use mcm_finder_lib::dataset::{Dataset, EndsCachedDataset, VecDataset};

fn criterion_benchmark(c: &mut Criterion) {
    {
        let mut c = c.benchmark_group("Basic VecDataset");
        c.sample_size(10).measurement_time(Duration::from_secs(30));

        c.bench_function("Loading MNIST11", |b| {
            b.iter(|| VecDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap())
        });
        c.bench_function("Loading Big5", |b| {
            b.iter(|| VecDataset::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap())
        });

        let dataset = VecDataset::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap();
        let icc = FixedBitSet::with_capacity_and_blocks(
            50,
            [0b0000000001000000000100000001000000000100000000010],
        );
        c.bench_function("Calculate Log E Big5", |b| b.iter(|| dataset.log_e(&icc)));
    }
    {
        let mut c = c.benchmark_group("Basic EndsCachedDataset");
        c.sample_size(10).measurement_time(Duration::from_secs(30));

        c.bench_function("Loading MNIST11", |b| {
            b.iter(|| {
                EndsCachedDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap()
            })
        });
        c.bench_function("Loading Big5", |b| {
            b.iter(|| {
                EndsCachedDataset::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap()
            })
        });

        let dataset = VecDataset::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap();
        let icc = FixedBitSet::with_capacity_and_blocks(
            50,
            [0b0000000001000000000100000001000000000100000000010],
        );
        c.bench_function("Calculate Log E Big5", |b| b.iter(|| dataset.log_e(&icc)));
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
            EndsCachedDataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap();
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
