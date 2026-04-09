use bit_vec::BitVec;
use statrs::function::gamma::ln_gamma;
use std::{f64::consts::PI, num::NonZeroU32};

/// Calculates the geometric complexity of an Independent Complete Component.
///
/// The geometric complexity of a Minimally Complex Model is equal to the
/// sum of the geometric complexity of its ICCs, with the complexity of an ICC being:
///
/// > 2ʳᵃ⁻¹ log π - log Γ(2ʳᵃ⁻¹)
///
/// with rₐ the modeled spin variables.
///
/// # Example
/// ```
/// # use mcm_finder::geometric_complexity_icc;
/// # use std::num::NonZeroU32;
/// let spin_variables = NonZeroU32::new(9).unwrap();
/// assert_eq!(geometric_complexity_icc(spin_variables), -868.6612503409542);
/// ```
pub fn geometric_complexity_icc(spin_variables: NonZeroU32) -> f64 {
    let spin_variables = u32::from(spin_variables);
    let pow: f64 = dbg!(2f64.powi((spin_variables - 1) as i32));
    dbg!(PI.ln() * pow) - dbg!(ln_gamma(pow))
}

/// Calculates the parameter complexity of an Independent Complete Component.
///
/// # Example
/// ```
/// # use mcm_finder::parameter_complexity_icc;
/// # use std::num::NonZeroU32;
/// let spin_variables = NonZeroU32::new(9).unwrap();
/// assert_eq!(parameter_complexity_icc(spin_variables, 5), -58.3662038406751);
/// ```
pub fn parameter_complexity_icc(spin_variables: NonZeroU32, n: u32) -> f64 {
    let spin_variables = u32::from(spin_variables);
    let K = 2f64.powi(spin_variables as i32) - 1.0;
    K * (((n as f64) / 2f64) / PI).ln() / 2.0
}

pub fn complexity_mcm(partition: Vec<BitVec>, n: u32, c_param: &mut f64, c_geom: &mut f64) -> f64 {
    *c_param = 0.0;
    *c_geom = 0.0;

    for part in partition.iter() {
        let spin_variables = NonZeroU32::try_from(part.count_ones() as u32).unwrap();
        *c_param += parameter_complexity_icc(spin_variables, n)
    }

    *c_param + *c_geom
}

/// The type for a Minimally Complex Model
///
/// # Example
/// ```
/// # use mcm_finder::MinimallyComplexModel;
/// # use bit_vec::BitVec;
/// let mcm = MinimallyComplexModel::new(
///     vec![
///         BitVec::from_bytes(&[0b11011100, 0b0]),
///         BitVec::from_bytes(&[0b00100011, 0b0]),
///         BitVec::from_bytes(&[0b00000000, 0b1]),
///     ],
///     9,
/// );
/// assert_eq!(mcm.complexity_mcm(), 1.3555732128424305);
/// ```
#[derive(Debug)]
pub struct MinimallyComplexModel {
    size: u32,
    basis: Vec<BitVec>,
}

impl MinimallyComplexModel {
    pub fn new(basis: Vec<BitVec>, size: u32) -> MinimallyComplexModel {
        MinimallyComplexModel { size, basis }
    }

    pub fn complexity_mcm(&self) -> f64 {
        let mut c_param = 0.0f64;
        let mut c_geom = 0.0f64;

        for part in self.basis.iter() {
            let spin_variables = NonZeroU32::try_from(part.count_ones() as u32).unwrap();
            c_param += parameter_complexity_icc(spin_variables, self.size);
            c_geom += geometric_complexity_icc(spin_variables);
        }

        c_param + c_geom
    }

    pub fn log_e(&self) -> f64 {}
}
