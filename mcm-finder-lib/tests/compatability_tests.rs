use std::num::NonZeroU32;

use approx::assert_relative_eq;
use mcm_finder_lib::{
    geometric_complexity_icc, parameter_complexity_icc,
};

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

// #[test]
// fn convert_to_mu_basis() {
//     let test_data: Vec<(usize, usize)> = vec![
//         (0b000000000, 0b000000000),
//         (0b000000001, 0b010000100),
//         (0b000000100, 0b100010001),
//         (0b000010000, 0b000001000),
//         (0b000010100, 0b100011001),
//         (0b000100000, 0b000000001),
//         (0b000100001, 0b010000101),
//         (0b000100100, 0b100010000),
//         (0b000101100, 0b100010010),
//         (0b000110100, 0b100011000),
//         (0b001000000, 0b010100000),
//         (0b001000001, 0b000100100),
//         (0b001000011, 0b001100010),
//         (0b001000100, 0b110110001),
//         (0b001001011, 0b001100000),
//         (0b001001100, 0b110110011),
//         (0b001010010, 0b011101110),
//         (0b001100100, 0b110110000),
//         (0b001101111, 0b101110000),
//         (0b001110100, 0b110111000),
//         (0b010000000, 0b001000000),
//         (0b010000001, 0b011000100),
//         (0b010000010, 0b000000110),
//         (0b010000100, 0b101010001),
//         (0b010000101, 0b111010101),
//         (0b010001011, 0b010000000),
//         (0b010001110, 0b100010101),
//         (0b010100000, 0b001000001),
//         (0b010100100, 0b101010000),
//         (0b010100110, 0b100010110),
//         (0b010101011, 0b010000001),
//         (0b010101100, 0b101010010),
//         (0b010101111, 0b110010000),
//         (0b010110000, 0b001001001),
//         (0b010110100, 0b101011000),
//         (0b011000011, 0b000100010),
//         (0b011001011, 0b000100000),
//         (0b011010011, 0b000101010),
//         (0b011010101, 0b101111101),
//         (0b011011001, 0b001101110),
//         (0b011011101, 0b101111111),
//         (0b011100111, 0b100110010),
//         (0b011101011, 0b000100001),
//         (0b011101101, 0b101110110),
//         (0b011111110, 0b110111100),
//         (0b100000000, 0b000111000),
//         (0b100000010, 0b001111110),
//         (0b100000100, 0b100101001),
//         (0b100010011, 0b011110010),
//         (0b100010100, 0b100100001),
//         (0b100011111, 0b111100001),
//         (0b100100000, 0b000111001),
//         (0b100100011, 0b011111011),
//         (0b100100100, 0b100101000),
//         (0b100100110, 0b101101110),
//         (0b100101100, 0b100101010),
//         (0b100110001, 0b010110101),
//         (0b100110100, 0b100100000),
//         (0b100110101, 0b110100100),
//         (0b100110110, 0b101100110),
//         (0b100110111, 0b111100010),
//         (0b100111100, 0b100100010),
//         (0b100111110, 0b101100100),
//         (0b101000000, 0b010011000),
//         (0b101000001, 0b000011100),
//         (0b101000011, 0b001011010),
//         (0b101000100, 0b110001001),
//         (0b101000101, 0b100001101),
//         (0b101001000, 0b010011010),
//         (0b101001011, 0b001011000),
//         (0b101001101, 0b100001111),
//         (0b101010000, 0b010010000),
//         (0b101010001, 0b000010100),
//         (0b101010010, 0b011010110),
//         (0b101010101, 0b100000101),
//         (0b101011010, 0b011010100),
//         (0b101100000, 0b010011001),
//         (0b101100001, 0b000011101),
//         (0b101100100, 0b110001000),
//         (0b101100101, 0b100001100),
//         (0b101100110, 0b111001110),
//         (0b101100111, 0b101001010),
//         (0b101101100, 0b110001010),
//         (0b101101101, 0b100001110),
//         (0b101101111, 0b101001000),
//         (0b101110000, 0b010010001),
//         (0b101110001, 0b000010101),
//         (0b101110100, 0b110000000),
//         (0b101110101, 0b100000100),
//         (0b101110110, 0b111000110),
//         (0b101110111, 0b101000010),
//         (0b101111011, 0b001010001),
//         (0b101111100, 0b110000010),
//         (0b101111101, 0b100000110),
//         (0b101111110, 0b111000100),
//         (0b101111111, 0b101000000),
//         (0b110000011, 0b010111010),
//         (0b110010011, 0b010110010),
//         (0b110011001, 0b011110110),
//         (0b110011011, 0b010110000),
//         (0b110011111, 0b110100001),
//         (0b110101110, 0b100101100),
//         (0b110110100, 0b101100000),
//         (0b110110110, 0b100100110),
//         (0b110110111, 0b110100010),
//         (0b110111110, 0b100100100),
//         (0b111000100, 0b111001001),
//         (0b111001011, 0b000011000),
//         (0b111010000, 0b011010000),
//         (0b111010011, 0b000010010),
//         (0b111010100, 0b111000001),
//         (0b111010101, 0b101000101),
//         (0b111010111, 0b100000011),
//         (0b111011000, 0b011010010),
//         (0b111011001, 0b001010110),
//         (0b111011011, 0b000010000),
//         (0b111011110, 0b110000101),
//         (0b111011111, 0b100000001),
//         (0b111100111, 0b100001010),
//         (0b111101101, 0b101001110),
//         (0b111110100, 0b111000000),
//         (0b111110101, 0b101000100),
//         (0b111110110, 0b110000110),
//         (0b111110111, 0b100000010),
//         (0b111111100, 0b111000010),
//         (0b111111101, 0b101000110),
//         (0b111111110, 0b110000100),
//         (0b111111111, 0b100000000),
//     ];
//     let mcm = MinimallyComplexModel::new(vec![
//         FixedBitSet::with_capacity_and_blocks(9, [0b000100100]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b000001010]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b000000011]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b100010000]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b100000100]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b101000000]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b010000010]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b001000001]),
//         FixedBitSet::with_capacity_and_blocks(9, [0b000000100]),
//     ]);
//     let dataset = Dataset::new(
//         test_data
//             .iter()
//             .map(|(i, _)| (FixedBitSet::with_capacity_and_blocks(9, [*i]), 1usize))
//             .collect(),
//         test_data.len(),
//     );
//     let correct_dataset = Dataset::new(
//         test_data
//             .iter()
//             .map(|(_, o)| (FixedBitSet::with_capacity_and_blocks(9, [*o]), 1usize))
//             .collect(),
//         test_data.len(),
//     );
//     let result_dataset = mcm.create_kset(dataset);
//     assert_eq!(result_dataset.data, correct_dataset.data);
// }

#[cxx::bridge]
mod ffi {
    unsafe extern "C++" {
        include!("mcm-finder-lib/cpp_functions/Complexity.cpp");

        fn GeomComplexity_ICC(m: u32) -> f64;

        fn ParamComplexity_ICC(m: u32, N: u32) -> f64;
    }
}
