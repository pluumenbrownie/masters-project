use std::{collections::HashMap, hash::BuildHasher, marker::PhantomData};

use fixedbitset::FixedBitSet;

pub trait DataContainer<S: BuildHasher + Default>: From<HashMap<FixedBitSet, usize, S>> {
    fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)>;

    fn name() -> String;

    fn to_icc(&self, icc: &FixedBitSet) -> Self {
        let mut new_icc_map = HashMap::with_hasher(S::default());
        for (o, n) in self.iter().map(|(i, n)| (i & icc, n)) {
            new_icc_map.entry(o).and_modify(|v| *v += n).or_insert(n);
        }
        Self::from(new_icc_map)
    }

    fn bins(&self) -> usize {
        self.iter().len()
    }
}

#[derive(Debug, Clone)]
pub struct DataVec<S: BuildHasher> {
    pub(crate) data: Vec<(FixedBitSet, usize)>,
    pub(crate) _hasher: PhantomData<S>,
}

impl<S: BuildHasher + Default> From<HashMap<FixedBitSet, usize, S>> for DataVec<S> {
    fn from(value: HashMap<FixedBitSet, usize, S>) -> Self {
        Self {
            data: value.into_iter().collect(),
            _hasher: PhantomData,
        }
    }
}

impl<S: BuildHasher + Default> DataContainer<S> for DataVec<S> {
    fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)> {
        self.data.iter().map(|(v, n)| (v, *n))
    }

    fn name() -> String {
        "DataVec".into()
    }
}

#[derive(Debug, Clone)]
pub struct DataMap<S: BuildHasher> {
    pub(crate) data: HashMap<FixedBitSet, usize, S>,
}

impl<S: BuildHasher> From<HashMap<FixedBitSet, usize, S>> for DataMap<S> {
    fn from(value: HashMap<FixedBitSet, usize, S>) -> Self {
        Self { data: value }
    }
}

impl<S: BuildHasher + Default> DataContainer<S> for DataMap<S> {
    fn iter(&self) -> impl ExactSizeIterator<Item = (&FixedBitSet, usize)> {
        self.data.iter().map(|(v, n)| (v, *n))
    }

    fn name() -> String {
        "DataMap".into()
    }
}
