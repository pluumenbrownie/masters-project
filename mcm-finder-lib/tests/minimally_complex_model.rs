use bit_vec::BitVec;
use mcm_finder_lib::MinimallyComplexModel;

#[test]
fn mcm() {
    let mcm = MinimallyComplexModel::new(
        vec![
            BitVec::from_bytes(&[0b11011100, 0b0]),
            BitVec::from_bytes(&[0b00100011, 0b0]),
            BitVec::from_bytes(&[0b00000000, 0b1]),
        ],
        9,
    );
    assert_eq!(mcm.complexity_mcm(), 1.3555732128424305);
}
