//! This module contains all of the solvers and relevant traits in this library.

mod annealing;
mod exhaustive;
mod greedy;
mod solvers_base;

pub use annealing::*;
pub use exhaustive::*;
pub use greedy::*;
pub use solvers_base::*;
