use miette::{Diagnostic, NamedSource, SourceSpan};
use thiserror::Error;

#[derive(Debug, Diagnostic, Error)]
pub enum MCMError {
    #[error("IO error: {0}")]
    #[diagnostic(code(mcm_finder_lib::io_error))]
    Io(#[from] std::io::Error),
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

    #[error("New basis overlaps with existing basis")]
    #[diagnostic(
        help("Elements should be in only one basis."),
        code("mcm-finder-lib::MCMError::BadLength")
    )]
    OverlappingBasis,
}
