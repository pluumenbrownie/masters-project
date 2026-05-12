use std::path::Path;

use miette::Result;

use mcm_finder_lib::solvers::{EvolutionarySolver, SelectionType, Solver};

fn main() -> Result<()> {
    // let filepath = Path::new("mcm-finder-lib/tests/data/SCOTUS_n9_N895_Data.dat");
    let filepath = Path::new("mcm-finder-lib/tests/data/MNIST11.sorted");
    // let filepath = Path::new("mcm-finder-lib/tests/data/Big5PT.sorted");

    let solver = EvolutionarySolver::from_file(filepath)?
        .set_generation_size(100)
        .set_shuffle_steps(50)
        .set_mutation_rate(1)
        .set_parent_selection(SelectionType::Exponential)
        .set_survivor_selection(SelectionType::Exponential);
    let result = solver.solve();
    println!("{result}");
    Ok(())
}
