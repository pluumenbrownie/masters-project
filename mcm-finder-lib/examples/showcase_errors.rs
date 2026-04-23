use std::path::Path;

use mcm_finder_lib::dataset::VecDataset;
use miette::{ErrReport, Result};

fn main() -> Result<()> {
    {
        println!("\nWhen data contains a bad character:");
        let dataset = VecDataset::read_from_file(Path::new(
            "mcm-finder-lib/tests/data/SCOTUS_n9_N11_bad_data.dat",
        ));
        if let Err(e) = dataset {
            println!("{:?}", ErrReport::from(e))
        };
    }
    {
        println!("\nWhen data contains a non-ASCII character:");
        let dataset = VecDataset::read_from_file(Path::new(
            "mcm-finder-lib/tests/data/SCOTUS_n9_N11_nonascii.dat",
        ));
        if let Err(e) = dataset {
            println!("{:?}", ErrReport::from(e))
        };
    }
    {
        println!("\nWhen data contains a row with the wrong length:");
        let dataset = VecDataset::read_from_file(Path::new(
            "mcm-finder-lib/tests/data/SCOTUS_n9_N11_bad_length.dat",
        ));
        if let Err(e) = dataset {
            println!("{:?}", ErrReport::from(e))
        };
    }
    Ok(())
}
