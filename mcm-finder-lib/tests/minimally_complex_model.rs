use fixedbitset::FixedBitSet;
use mcm_finder_lib::MinimallyComplexModel;

#[test]
fn mcm() {
    let mcm = MinimallyComplexModel::new(
        vec![
            FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
            FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
            FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
        ],
        9,
    );
    assert_eq!(mcm.complexity_mcm(), 1.3555732128424305);
}
