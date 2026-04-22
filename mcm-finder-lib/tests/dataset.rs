use std::path::Path;

use fixedbitset::FixedBitSet;
use mcm_finder_lib::dataset::Dataset;
use pretty_assertions::assert_eq;

#[test]
fn partition() {
    let dataset =
        Dataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat")).unwrap();

    let icc = FixedBitSet::with_capacity_and_blocks(9, [0b001011101]);

    let partitioned_dataset = dataset.partition(&icc);

    assert_eq!(partitioned_dataset.datapoints(), dataset.datapoints());
    assert_eq!(
        partitioned_dataset.datapoints(),
        partitioned_dataset.iter().map(|p| p.1).sum()
    );
    assert_eq!(partitioned_dataset.get(&icc), Some(441));
    let should_fit_in_icc = vec![
        FixedBitSet::with_capacity_and_blocks(9, [0b001011101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b001011111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b001111101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b001111111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b011011101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b011011111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b011111101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b011111111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b101011101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b101011111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b101111101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b101111111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b111011101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b111011111]),
        FixedBitSet::with_capacity_and_blocks(9, [0b111111101]),
        FixedBitSet::with_capacity_and_blocks(9, [0b111111111]),
    ];
    let test_total: usize = should_fit_in_icc
        .iter()
        .map(|v| dataset.get(v).unwrap_or(0))
        .sum();
    assert_eq!(Some(test_total), partitioned_dataset.get(&icc))
}
