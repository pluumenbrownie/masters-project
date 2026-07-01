use fixedbitset::FixedBitSet;

/// The range of an ICC which contains included variables.
///
/// - `start` is the index of the first included variable.
///
/// - `end` is the index of the first excluded variable after the last included variable.
#[derive(Debug, Clone)]
pub(crate) struct Ends {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl Ends {
    pub(crate) fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    pub(crate) fn from_icc(partition: &FixedBitSet) -> Option<Self> {
        partition.minimum().map(|start| Ends {
            start,
            end: partition.maximum().unwrap() + 1,
        })
    }
}

pub(crate) const fn get_ends_cache_location(ends: Option<Ends>, variables: usize) -> usize {
    match ends {
        None => (variables * (variables + 1)) / 2,
        Some(ends) => {
            if ends.start == 0 && ends.end == variables + 1 {
                0
            } else {
                let start_offset = variables - ends.start;
                (variables.pow(2) + variables - (start_offset.pow(2) + start_offset)) / 2
                    + variables
                    - ends.end
            }
        }
    }
}
