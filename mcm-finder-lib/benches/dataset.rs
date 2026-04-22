use std::path::Path;

use criterion::{Criterion, criterion_group, criterion_main};

use fixedbitset::FixedBitSet;
use mcm_finder_lib::dataset::Dataset;

fn criterion_benchmark(c: &mut Criterion) {
    {
        let dataset = Dataset::read_from_file(Path::new("tests/data/MNIST11.sorted")).unwrap();
        let icc = FixedBitSet::with_capacity_and_blocks(
            122,
            [
                0b0000000001000000000100000001000000000100000000010,
                0b000000000100000000010000000100000000,
            ],
        );
        c.bench_function("Dataset partition MNIST", |b| {
            b.iter(|| dataset.partition(&icc))
        });
    }
    {
        let dataset = Dataset::read_from_file(Path::new("tests/data/Big5PT.sorted")).unwrap();
        let icc = FixedBitSet::with_capacity_and_blocks(
            50,
            [0b0000000001000000000100000001000000000100000000010],
        );
        c.bench_function("Dataset partition Big5", |b| {
            b.iter(|| dataset.partition(&icc))
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
