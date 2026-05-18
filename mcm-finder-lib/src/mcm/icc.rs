use std::{
    fmt::Display,
    ops::{BitAnd, BitAndAssign, BitOrAssign, Deref, DerefMut},
    sync::OnceLock,
};

use fixedbitset::FixedBitSet;

use crate::{
    dataset::{Dataset, LogE},
    mcm::gamma_factor,
};

#[derive(Debug)]
pub struct IndependentCompleteComponent {
    pub(crate) bits: FixedBitSet,
    pub(crate) log_e: OnceLock<f64>,
}

impl IndependentCompleteComponent {
    pub fn log_e<T: Dataset>(&self, dataset: &T) -> f64 {
        *self.log_e.get_or_init(|| self.calculate_log_e(dataset))
    }

    pub(crate) fn calculate_log_e<T: Dataset>(&self, dataset: &T) -> f64 {
        let rank_subset: i32 = self.count_ones(..).try_into().unwrap();

        let gamma_factor = gamma_factor(dataset, rank_subset);

        gamma_factor + dataset.transform_to_icc(self).log_e()
    }

    pub fn full_clone(&self) -> Self {
        IndependentCompleteComponent {
            bits: self.bits.clone(),
            log_e: self.log_e.clone(),
        }
    }
}

impl Clone for IndependentCompleteComponent {
    fn clone(&self) -> Self {
        IndependentCompleteComponent {
            bits: self.bits.clone(),
            log_e: OnceLock::new(),
        }
    }
}

impl Display for IndependentCompleteComponent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.bits)
    }
}

impl PartialEq for IndependentCompleteComponent {
    fn eq(&self, other: &Self) -> bool {
        self.bits == other.bits
    }
}

impl Eq for IndependentCompleteComponent {}

impl From<FixedBitSet> for IndependentCompleteComponent {
    fn from(value: FixedBitSet) -> Self {
        Self {
            bits: value,
            log_e: OnceLock::new(),
        }
    }
}

impl Deref for IndependentCompleteComponent {
    type Target = FixedBitSet;
    fn deref(&self) -> &Self::Target {
        &self.bits
    }
}

impl DerefMut for IndependentCompleteComponent {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.log_e = OnceLock::new();
        &mut self.bits
    }
}

impl BitAnd for &IndependentCompleteComponent {
    type Output = IndependentCompleteComponent;
    fn bitand(self, rhs: Self) -> Self::Output {
        IndependentCompleteComponent::from(&self.bits & &rhs.bits)
    }
}

impl BitOrAssign<&Self> for IndependentCompleteComponent {
    fn bitor_assign(&mut self, rhs: &Self) {
        self.log_e = OnceLock::new();
        self.bits |= &rhs.bits;
    }
}

impl BitAndAssign<&Self> for IndependentCompleteComponent {
    fn bitand_assign(&mut self, rhs: &Self) {
        self.log_e = OnceLock::new();
        self.bits &= &rhs.bits;
    }
}

impl BitAndAssign<Self> for IndependentCompleteComponent {
    fn bitand_assign(&mut self, rhs: Self) {
        self.log_e = OnceLock::new();
        self.bits &= &rhs.bits;
    }
}
