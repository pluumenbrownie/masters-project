use std::path::Path;

use approx::assert_relative_eq;
use fixedbitset::FixedBitSet;
use mcm_finder_lib::dataset::{
    DataContainter, DataMap, DataVec, Dataset, EndsCachedDataset, SimpleDataset,
};
use pretty_assertions::assert_eq;

#[test]
fn ends_cache_partition_vec() {
    ends_cache_partition::<DataVec>();
}

#[test]
fn ends_cache_partition_map() {
    ends_cache_partition::<DataMap>();
}

fn ends_cache_partition<C: DataContainter>() {
    let ends_dataset =
        EndsCachedDataset::<C>::read_from_file(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat"))
            .unwrap();
    let reference_dataset =
        SimpleDataset::<C>::read_from_file(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat"))
            .unwrap();
    assert_eq!(
        ends_dataset.data.len(),
        (ends_dataset.variables() * (ends_dataset.variables() + 1)) / 2 + 1
    );

    let icc = FixedBitSet::with_capacity_and_blocks(9, [0b001011101]);
    let partitioned_dataset = ends_dataset.transform_to_icc(&icc);

    assert_eq!(partitioned_dataset.datapoints(), ends_dataset.datapoints());
    assert_eq!(
        partitioned_dataset.datapoints(),
        partitioned_dataset.iter().map(|p| p.1).sum()
    );
    println!("{partitioned_dataset}");
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
        .map(|v| ends_dataset.get_icc(v).unwrap_or(0))
        .sum();
    assert_eq!(Some(test_total), partitioned_dataset.get(&icc));

    assert_relative_eq!(
        partitioned_dataset.log_e(&icc),
        reference_dataset.log_e(&icc),
        epsilon = 1e-9,
    );

    assert_relative_eq!(
        ends_dataset.log_e(&FixedBitSet::with_capacity_and_blocks(9, [0])),
        5188.503754891344,
        epsilon = 1e-9,
    );
    assert_relative_eq!(
        ends_dataset.log_e(&FixedBitSet::with_capacity_and_blocks(9, [0b100000001])),
        4136.3735176022055,
        epsilon = 1e-9,
    );
}
