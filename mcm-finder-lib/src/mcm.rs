use fixedbitset::FixedBitSet;
use miette::NamedSource;
use statrs::function::gamma::ln_gamma;
use std::{
    collections::HashMap,
    f64::consts::{LN_2, PI},
    fmt::Display,
    fs::File,
    io::{BufReader, Read},
    num::{NonZeroU32, NonZeroUsize},
    path::Path,
};

use crate::{
    dataset::{Dataset, LogE, line_length_tracker, verify_ascii},
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
/// let mcm = MinimallyComplexModel::from_iccs(
///     vec![
///         FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
///         FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
///         FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
///     ],
/// ).unwrap();
/// assert_eq!(mcm.rank(), 9);
/// assert_eq!(mcm.complexity_mcm(), 1.3555732128424305);
/// ```
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct MinimallyComplexModel {
    pub(crate) partition: Vec<FixedBitSet>,
}

impl MinimallyComplexModel {
    /// Create a new MCM with sorted ICCs.
    fn new(partition: Vec<FixedBitSet>) -> MinimallyComplexModel {
        let mut partition = partition;
        partition.sort();
        MinimallyComplexModel { partition }
    }

    /// Returns `true` if ICCs in array do not overlap anywhere. The ICCs do not
    /// have to include all of the variables.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mut partition = vec!(
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// );
    /// assert!(MinimallyComplexModel::verify_iccs(&partition));
    /// partition[1].set(0, true);
    /// assert!(!MinimallyComplexModel::verify_iccs(&partition));
    /// ```
    pub fn verify_iccs(partition: &[FixedBitSet]) -> bool {
        for (nr, one) in partition.iter().enumerate() {
            for two in &partition[(nr + 1)..] {
                if !(one & two).is_clear() {
                    return false;
                }
            }
        }
        true
    }

    /// Returns the trivial MCM with each variable in a different partition. Note that
    /// ICCs in an MCM will be sorted so may have a different order from how they were
    /// inserted.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mut partition = vec!(
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// );
    /// let mcm = MinimallyComplexModel::from_iccs(partition).unwrap();
    /// assert_eq!(mcm.rank(), 9);
    /// assert_eq!(mcm.count_icc(), 3);
    /// ```
    pub fn from_iccs(partition: Vec<FixedBitSet>) -> Result<MinimallyComplexModel, MCMError> {
        if MinimallyComplexModel::verify_iccs(&partition) {
            let mut partition = partition;
            partition.sort();
            return Ok(MinimallyComplexModel { partition });
        }
        Err(MCMError::FromIccs)
    }

    /// Returns an MCM with a single partition.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::full(NonZero::new(100).unwrap());
    /// assert_eq!(mcm.rank(), 100);
    /// assert_eq!(mcm.count_icc(), 1);
    /// ```
    pub fn full(variables: NonZeroUsize) -> MinimallyComplexModel {
        let mut partition = FixedBitSet::with_capacity(variables.into());
        partition.set_range(.., true);
        MinimallyComplexModel {
            partition: vec![partition],
        }
    }

    /// Returns the trivial MCM with each variable in a different partition.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::trivial(NonZero::new(500).unwrap());
    /// assert_eq!(mcm.rank(), 500);
    /// assert_eq!(mcm.count_icc(), 500);
    /// ```
    pub fn trivial(variables: NonZeroUsize) -> MinimallyComplexModel {
        let variables: usize = variables.into();
        let mut partition = Vec::with_capacity(variables);
        let vector_length = variables.div_ceil(usize::BITS as usize);
        let mut part_content = vec![0usize; vector_length];
        let mut counter = 0usize;

        for n in 0..vector_length {
            part_content[n] = 1;
            for _ in 0..usize::BITS {
                partition.push(FixedBitSet::with_capacity_and_blocks(
                    variables,
                    part_content.clone(),
                ));
                part_content[n] <<= 1;

                counter += 1;
                if counter == variables {
                    break;
                }
            }
        }
        MinimallyComplexModel { partition }
    }

    /// Returns the rank of this MCM. Can be lower than `self.variables()` if not
    /// all variables are included in the model.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]).unwrap();
    /// assert_eq!(mcm.rank(), 9);
    /// ```
    /// Now if we remove the first variable from this model:
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use fixedbitset::FixedBitSet;
    /// let smaller_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b010111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]).unwrap();
    /// assert_eq!(smaller_mcm.rank(), 8);
    /// ```
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
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    /// ]).unwrap();
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
    /// let mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b000000001]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111000]),
    /// ]).unwrap();
    /// let result_mcm = MinimallyComplexModel::from_iccs(vec![
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b001000110]),
    ///     FixedBitSet::with_capacity_and_blocks(9, [0b110111001]),
    /// ]).unwrap();
    /// assert_eq!(mcm.merge(2, 0), result_mcm);
    /// ```
    pub fn merge(&self, basis: usize, into: usize) -> MinimallyComplexModel {
        let mut partition = self.partition.clone();

        partition[into] |= &self.partition[basis];
        partition.remove(basis);

        MinimallyComplexModel::new(partition)
    }

    /// Returns the amount of ICCs present in this model.
    ///
    /// # Example
    /// ```
    /// # use mcm_finder_lib::mcm::MinimallyComplexModel;
    /// # use std::num::NonZero;
    /// let mcm = MinimallyComplexModel::trivial(NonZero::new(5).unwrap());
    /// assert_eq!(mcm.count_icc(), 5);
    /// ```
    pub fn count_icc(&self) -> usize {
        self.partition.len()
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
    pub fn log_e<T: Dataset>(
        &self,
        dataset: &T,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
    ) -> f64 {
        let mut log_e = 0f64;

        for part in self.partition.iter() {
            let rank_subset: i32 = part.count_ones(..).try_into().unwrap();

            let gamma_factor = ln_gamma(2.0f64.powi(rank_subset - 1))
                - ln_gamma(dataset.datapoints() as f64 + 2.0f64.powi(rank_subset - 1));

            let sum_of_partitions = if let Some(cache) = log_e_cache {
                *cache
                    .entry(part.clone())
                    .or_insert_with(|| dataset.transform_to_icc(part).log_e())
            } else {
                dataset.transform_to_icc(part).log_e()
            };

            log_e += gamma_factor;
            log_e += sum_of_partitions;
        }

        let front_constant: f64 =
            (dataset.datapoints() * (self.variables() - self.rank())) as f64 * LN_2;
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

        Ok(MinimallyComplexModel { partition: vectors })
    }
}

impl Display for MinimallyComplexModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut partition_string: String = self
            .partition
            .iter()
            .flat_map(|p| format!("{p}\n").chars().collect::<Vec<_>>())
            .collect();
        partition_string.pop();
        write!(f, "{}", partition_string)
    }
}
