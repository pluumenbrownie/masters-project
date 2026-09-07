use std::path::Path;

use fixedbitset::FixedBitSet;

use mcm_finder_lib::dataset::{Dataset, simple::VecDataset};

#[test]
fn read_data() {
    let dataset = VecDataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat"));
    // println!("{:?}", dataset);
    assert!(dataset.is_ok());
    let unwrapped_dataset = dataset.unwrap();
    for (line, nr) in unwrapped_dataset.iter() {
        println!("{line} - {nr}");
    }

    assert_eq!(
        unwrapped_dataset.get(&FixedBitSet::with_capacity_and_blocks(9, vec![0b111111111])),
        Some(221)
    );
    assert_eq!(
        unwrapped_dataset.get(&FixedBitSet::with_capacity_and_blocks(9, vec![0])),
        Some(174)
    );
    assert_eq!(unwrapped_dataset.variable_states(), 2);
    assert_eq!(unwrapped_dataset.variable_width(), 1);
}

#[test]
fn read_empty() {
    let dataset = VecDataset::read_from_file(Path::new("./tests/data/empty.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_err())
}

#[test]
fn read_bad_length() {
    let dataset =
        VecDataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N11_bad_length.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_err())
}

#[test]
fn read_nonascii() {
    let dataset = VecDataset::read_from_file(Path::new("./tests/data/SCOTUS_n9_N11_nonascii.dat"));
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
