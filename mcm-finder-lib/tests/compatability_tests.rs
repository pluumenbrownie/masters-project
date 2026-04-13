use std::num::NonZeroU32;

use approx::assert_relative_eq;
use mcm_finder_lib::{geometric_complexity_icc, parameter_complexity_icc};
use quickcheck_macros::quickcheck;

use crate::ffi::{GeomComplexity_ICC, ParamComplexity_ICC};

#[test]
/// The C++ version of lgamma is not accurate beyond `spin_variables=64`.
fn geo_compl() {
    for spin_variables in 1u32..=64 {
        let epsilon: f64 = (GeomComplexity_ICC(spin_variables) / 1e9).abs();
        assert_relative_eq!(
            geometric_complexity_icc(NonZeroU32::try_from(spin_variables).unwrap()),
            GeomComplexity_ICC(spin_variables),
            epsilon = epsilon
        );
    }
}

#[test]
/// The C++ version is not accurate beyond `spin_variables=32`.
fn par_compl() {
    for spin_variables in 1u32..=32 {
        let epsilon: f64 = (GeomComplexity_ICC(spin_variables) / 1e9).abs();

        let spin_variables = NonZeroU32::try_from(spin_variables).unwrap();
        println!("{spin_variables}");

        assert_relative_eq!(
            dbg!(parameter_complexity_icc(spin_variables, 5)),
            dbg!(ParamComplexity_ICC(spin_variables.into(), 5)),
            epsilon = epsilon
        );
    }
}

#[cxx::bridge]
mod ffi {
    unsafe extern "C++" {
        include!("mcm-finder-lib/cpp_functions/Complexity.cpp");

        fn GeomComplexity_ICC(m: u32) -> f64;

        fn ParamComplexity_ICC(m: u32, N: u32) -> f64;
    }
}
