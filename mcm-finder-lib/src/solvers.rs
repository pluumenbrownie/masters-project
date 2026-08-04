//! This module contains all of the solvers and relevant traits in this library.

mod annealing;
mod ant_colony;
mod constructive;
mod custom_annealing;
mod evolutionary;
mod exhaustive;
mod greedy;
mod hill_climber;
mod parallel_tempering;
mod pso;
mod solvers_base;
mod spso;

pub use annealing::anneal_temps::*;
pub use annealing::*;
pub use ant_colony::*;
pub use constructive::*;
pub use custom_annealing::*;
pub use evolutionary::*;
pub use exhaustive::*;
pub use greedy::*;
pub use hill_climber::*;
pub use parallel_tempering::*;
pub use pso::*;
pub use solvers_base::*;
pub use spso::*;
