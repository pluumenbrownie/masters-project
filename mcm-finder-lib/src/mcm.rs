use fixedbitset::FixedBitSet;
use miette::NamedSource;
use statrs::function::gamma::ln_gamma;
use std::{
    f64::consts::{LN_2, PI},
    fmt::Display,
    fs::File,
    io::{BufReader, Read},
    num::NonZeroU32,
    path::Path,
};

use crate::{
    dataset::{Dataset, line_length_tracker, verify_ascii},
    mcm_error::MCMError,
};

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
/// # use mcm_finder_lib::mcm::geometric_complexity_icc;
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
/// # use mcm_finder_lib::mcm::parameter_complexity_icc;
/// # use std::num::NonZeroU32;
/// let spin_variables = NonZeroU32::new(9).unwrap();
/// assert_eq!(parameter_complexity_icc(spin_variables, 5), -58.3662038406751);
/// ```
pub fn parameter_complexity_icc(spin_variables: NonZeroU32, n: usize) -> f64 {
    let spin_variables = u32::from(spin_variables);
    let k = 2f64.powi(spin_variables as i32) - 1.0;
    k * (((n as f64) / 2f64) / PI).ln() / 2.0
}

/// The type for a Minimally Complex Model
///
/// # Example
/// ```
/// # use mcm_finder_lib::mcm::MinimallyComplexModel;
/// # use fixedbitset::FixedBitSet;
/// let mcm = MinimallyComplexModel::new(
///     vec![
///         FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
///         FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
///         FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
///     ],
/// );
/// assert_eq!(mcm.rank(), 9);
/// assert_eq!(mcm.complexity_mcm(), 1.3555732128424305);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MinimallyComplexModel {
    pub(crate) partition: Vec<FixedBitSet>,
}

impl MinimallyComplexModel {
    pub fn new(partition: Vec<FixedBitSet>) -> MinimallyComplexModel {
        MinimallyComplexModel { partition }
    }

    /// Returns the rank of this MCM. Can be lower than `self.variables()` if not
    /// all variables are included in the model.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::new(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]);
    /// assert_eq!(mcm.rank(), 9);
    /// ```
    /// Now if we remove the first variable from this model:
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let smaller_mcm = MinimallyComplexModel::new(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b010111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]);
    /// assert_eq!(smaller_mcm.rank(), 8);```
    pub fn rank(&self) -> usize {
        self.partition
            .iter()
            .fold(FixedBitSet::with_capacity(self.variables()), |a, b| &a | b)
            .count_ones(..)
    }

    /// Returns the amount of variables in this MCM.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::new(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]);
    /// assert_eq!(mcm.variables(), 9);
    /// ```
    pub fn variables(&self) -> usize {
        self.partition[0].len()
    }

    /// Merge the first basis into the second bases, and return a new MCM.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::new(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]);
    /// let result_mcm = MinimallyComplexModel::new(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    /// ]);
    /// assert_eq!(mcm.merge(2, 0), result_mcm);
    /// ```
    pub fn merge(&self, basis: usize, into: usize) -> MinimallyComplexModel {
        let mut partition = self.partition.clone();

        partition[into] |= &self.partition[basis];
        partition.remove(basis);

        MinimallyComplexModel { partition }
    }

    pub fn complexity_mcm(&self) -> f64 {
        let mut c_param = 0.0f64;
        let mut c_geom = 0.0f64;

        for part in self.partition.iter() {
            let spin_variables = NonZeroU32::try_from(part.count_ones(..) as u32).unwrap();
            c_param += parameter_complexity_icc(spin_variables, self.rank());
            c_geom += geometric_complexity_icc(spin_variables);
        }

        c_param + c_geom
    }

    /// Calculate the logarithm of the evidence of this MCM, via the equation
    ///
    /// log⁡E=-(n-r)ln2 + ∑_{P∈P}\[(ln⁡Γ(2ᴾ⁻¹)−ln⁡Γ(n+2ᴾ⁻¹))+∑_{x∈data}(P)(ln⁡Γ(kx+0.5)−ln⁡Γ(0.5))\]
    ///
    /// Hope that clears things up.
    pub fn log_e(&self, dataset: &Dataset) -> f64 {
        let mut log_e = 0f64;

        for part in self.partition.iter() {
            let rank_subset: i32 = part.count_ones(..).try_into().unwrap();

            let gamma_factor = ln_gamma(2.0f64.powi(rank_subset - 1))
                - ln_gamma(dataset.datapoints as f64 + 2.0f64.powi(rank_subset - 1));

            let sum_of_partitions = dataset
                .partition(part)
                .iter()
                .map(|(_, &k)| ln_gamma(k as f64 + 0.5) - ln_gamma(0.5))
                .sum::<f64>();

            log_e += gamma_factor;
            log_e += sum_of_partitions;
        }

        let front_constant: f64 =
            (dataset.datapoints * (self.variables() - self.rank())) as f64 * LN_2;
        log_e - front_constant
    }

    pub fn read_from_file(path: &Path) -> Result<MinimallyComplexModel, MCMError> {
        let mut vectors: Vec<FixedBitSet> = vec![];
        let filename = path.file_name().unwrap().to_str().unwrap().to_owned();
        let mut buf_reader = BufReader::new(File::open(path)?);

        // reading the entire file to a string
        let file = {
            let mut file = String::new();
            buf_reader.read_to_string(&mut file)?;
            file
        };

        let mut line_length = 0usize;
        let mut bool_array: Vec<bool> = vec![];

        for (nr, byte) in file.bytes().enumerate() {
            // validate character is valid ascii
            verify_ascii(&filename, &file, nr, byte)?;
            // count ones and zeroes
            match byte {
                b'0' => bool_array.push(false),
                b'1' => bool_array.push(true),
                b'\r' | b'\n' if !bool_array.is_empty() => {
                    // check the line length
                    line_length_tracker(&filename, &file, &mut line_length, &bool_array, nr)?;

                    // add the datapoints to the vector
                    let mut bitvec = FixedBitSet::with_capacity(bool_array.len());
                    for (nr, bit) in bool_array.drain(..).enumerate() {
                        bitvec.set(nr, bit);
                    }
                    vectors.push(bitvec);
                    debug_assert!(bool_array.is_empty());
                }
                b'\r' | b'\n' => {}
                // wrong character case
                _ => Err(MCMError::BadCharacter {
                    src: NamedSource::new(&filename, file.clone()),
                    bad_line: nr.into(),
                })?,
            }
        }

        Ok(MinimallyComplexModel::new(vectors))
    }
}

impl Display for MinimallyComplexModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let partition_string: String = self
            .partition
            .iter()
            .flat_map(|p| format!("\n{p}").chars().collect::<Vec<_>>())
            .collect();
        write!(f, "MCM with vectors: {}", partition_string)
    }
}
