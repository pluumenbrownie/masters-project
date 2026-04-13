use std::path::Path;

use mcm_finder_lib::Dataset;
use miette::Result;

fn main() -> Result<()> {
    let dataset = Dataset::read(Path::new(
        "mcm-finder-lib/tests/data/SCOTUS_n9_N11_bad_length.dat",
    ));
    println!("{:?}", dataset);
    dataset?;
    Ok(())
}
