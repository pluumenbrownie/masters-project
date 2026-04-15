use std::path::Path;

use fixedbitset::FixedBitSet;
use mcm_finder_lib::Dataset;

#[test]
fn read_data() {
    let dataset = Dataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_ok());

    let unwrapped_dataset = dataset.unwrap();

    assert_eq!(
        unwrapped_dataset.get(FixedBitSet::with_capacity_and_blocks(9, vec![0b111111111])),
        Some(221)
    );
    assert_eq!(
        unwrapped_dataset.get(FixedBitSet::with_capacity_and_blocks(9, vec![0])),
        Some(174)
    );
}

#[test]
fn read_bad_data() {
    let dataset = Dataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N11_bad_data.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_err())
}

#[test]
fn read_bad_length() {
    let dataset = Dataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N11_bad_length.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_err())
}

#[test]
fn read_nonascii() {
    let dataset = Dataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N11_nonascii.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_err())
}

// #[test]
// fn read_basis_vectors() {
//     let dataset = Dataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat"));
//     println!("{:?}", dataset);
//     assert!(dataset.is_ok());

//     let unwrapped_dataset = dataset.unwrap();

//     assert_eq!(
//         unwrapped_dataset.get(FixedBitSet::with_capacity_and_blocks(9, vec![0b111111111])),
//         Some(221)
//     );
//     assert_eq!(
//         unwrapped_dataset.get(FixedBitSet::with_capacity_and_blocks(9, vec![0])),
//         Some(174)
//     );
// }
