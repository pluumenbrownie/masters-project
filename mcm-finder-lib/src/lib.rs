use bit_vec::BitVec;
use miette::{Diagnostic, NamedSource, SourceSpan};
use statrs::function::gamma::ln_gamma;
use std::{
    collections::HashMap,
    f64::consts::PI,
    fs::File,
    io::{BufReader, Read},
    num::NonZeroU32,
    path::Path,
};
use thiserror::Error;

#[derive(Debug, Diagnostic, Error)]
pub enum MCMError {
    #[error("IO error: {0}")]
    #[diagnostic(code(mcm_finder_lib::io_error))]
    Io(#[from] std::io::Error),
    // #[error("Line has incorrect length")]
    // WrongLength {
    //     #[label("here")]
    //     line: NamedSource<String>,
    // },
    #[error("Bad character error")]
    #[diagnostic(
        help("Data lines should only consist of 0 and 1, and should be of equal length."),
        code("mcm-finder-lib::MCMError::BadCharacter")
    )]
    BadCharacter {
        #[source_code]
        src: NamedSource<String>,
        #[label("Incorrect character")]
        bad_line: SourceSpan,
    },
    #[error("Bad length error")]
    #[diagnostic(
        help("Data lines should have the same length."),
        code("mcm-finder-lib::MCMError::BadLength")
    )]
    BadLength {
        #[source_code]
        src: NamedSource<String>,
        #[label("Bad line")]
        bad_line: SourceSpan,
    },

    #[error("{filename} is empty")]
    #[diagnostic(help(
        "Data lines should only consist of 0 and 1, and should be of equal length."
    ))]
    EmptyFile { filename: String },
}

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
/// # use mcm_finder_lib::geometric_complexity_icc;
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
/// # use mcm_finder_lib::parameter_complexity_icc;
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
/// # use mcm_finder_lib::MinimallyComplexModel;
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

    // pub fn log_e(&self) -> f64 {}
}

#[derive(Debug)]
pub struct Dataset {
    data: HashMap<BitVec<usize>, usize>,
}

impl Dataset {
    pub fn new(data: HashMap<BitVec<usize>, usize>) -> Dataset {
        Dataset { data }
    }

    pub fn read(path: &Path) -> Result<Dataset, MCMError> {
        let mut data = HashMap::new();
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
            if !byte.is_ascii() {
                Err(MCMError::BadCharacter {
                    src: NamedSource::new(&filename, file.clone()),
                    bad_line: nr.into(),
                })?
            };
            // count ones and zeroes
            match byte {
                b'0' => bool_array.push(false),
                b'1' => bool_array.push(true),
                b'\r' | b'\n' if !bool_array.is_empty() => {
                    // check the line length
                    if line_length == 0 {
                        line_length = bool_array.len();
                    } else if bool_array.len() != line_length {
                        Err(MCMError::BadLength {
                            src: NamedSource::new(&filename, file.clone()),
                            bad_line: (nr - bool_array.len(), bool_array.len()).into(),
                        })?
                    }

                    // add the datapoints to the hashmap
                    let bitvec = BitVec::<usize>::from_iter(bool_array.drain(..));
                    data.entry(bitvec).and_modify(|i| *i += 1).or_insert(1usize);
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

        Ok(Dataset::new(data))
    }

    pub fn get(&self, bitvector: BitVec<usize>) -> Option<usize> {
        self.data.get(&bitvector).copied()
    }
}
