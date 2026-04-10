use std::path::Path;

use mcm_finder::Dataset;

#[test]
fn read_data() {
    let dataset = Dataset::read(Path::new("./tests/data/SCOTUS_n9_N895_Data.dat"));
    println!("{:?}", dataset);
    assert!(dataset.is_ok())
}

#[test]
fn read_bad_data() {
    {
        let dataset = Dataset::read(Path::new("./tests/data/SCOTUS_n9_N11_bad_data.dat"));
        println!("{:?}", dataset);
        assert!(dataset.is_err())
    }
    {
        let dataset = Dataset::read(Path::new("./tests/data/SCOTUS_n9_N11_bad_length.dat"));
        println!("{:?}", dataset);
        assert!(dataset.is_err())
    }
}
