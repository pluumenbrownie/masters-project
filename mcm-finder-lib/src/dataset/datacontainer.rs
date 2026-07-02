use std::{collections::HashMap, hash::BuildHasher, marker::PhantomData};

use fixedbitset::FixedBitSet;

pub trait DataContainer<S: BuildHasher + Default>: From<HashMap<FixedBitSet, usize, S>> {
    fn iter(&self) -> impl ExactSizeIterator<Item = (FixedBitSet, usize)>;

    fn name() -> String;

    fn to_icc(&self, icc: &FixedBitSet) -> Self {
        let mut new_icc_map = HashMap::with_hasher(S::default());
        for (o, n) in self.iter().map(|(mut i, n)| {
            i &= icc;
            (i, n)
        }) {
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
    fn iter(&self) -> impl ExactSizeIterator<Item = (FixedBitSet, usize)> {
        self.data.iter().map(|(v, n)| (v.clone(), *n))
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
    fn iter(&self) -> impl ExactSizeIterator<Item = (FixedBitSet, usize)> {
        self.data.iter().map(|(v, n)| (v.clone(), *n))
    }

    fn name() -> String {
        "DataMap".into()
    }
}

#[derive(Debug, Default)]
pub struct DataVecSoA<S: BuildHasher> {
    data: Vec<Vec<usize>>,
    variables: usize,
    counts: Vec<usize>,
    _hasher: PhantomData<S>,
}

impl<S: BuildHasher> DataVecSoA<S> {
    fn init(&mut self, variables: usize) {
        assert_eq!(self.variables, 0);
        let bins = vec![Vec::default(); variables.div_ceil(usize::BITS as usize)];
        self.data = bins;
        self.variables = variables;
        self.counts = Vec::default();
    }

    pub fn push(&mut self, value: (FixedBitSet, usize)) {
        if self.variables == 0 {
            self.init(value.0.len());
        }
        assert_eq!(self.variables, value.0.len());
        for (v, x) in self.data.iter_mut().zip(value.0.as_slice()) {
            v.push(*x);
        }
        self.counts.push(value.1);
    }

    fn mask(&self, icc: &FixedBitSet) -> DataVecSoA<S> {
        let new_data = self
            .data
            .iter()
            .zip(icc.as_slice())
            .map(|(d, i)| d.iter().map(|x| x & i).collect())
            .collect();
        Self {
            data: new_data,
            variables: self.variables,
            counts: self.counts.clone(),
            _hasher: self._hasher,
        }
    }
}

impl<S: BuildHasher + Default> From<(FixedBitSet, usize)> for DataVecSoA<S> {
    fn from(value: (FixedBitSet, usize)) -> Self {
        let variables = value.0.len();
        let mut output = Self::default();
        output.init(variables);
        output.push(value);
        output
    }
}

impl<S: BuildHasher + Default> From<HashMap<FixedBitSet, usize, S>> for DataVecSoA<S> {
    fn from(value: HashMap<FixedBitSet, usize, S>) -> Self {
        let mut output = Self::default();
        for t in value.into_iter() {
            output.push(t);
        }
        output
    }
}

impl<S: BuildHasher + Default> DataContainer<S> for DataVecSoA<S> {
    fn bins(&self) -> usize {
        if self.variables == 0 {
            return 0;
        };
        self.data[0].len()
    }

    fn name() -> String {
        "DataVecSoA".into()
    }

    fn iter(&self) -> impl ExactSizeIterator<Item = (FixedBitSet, usize)> {
        (0..self.data[0].len())
            .zip(self.counts.iter())
            .map(|(i, n)| {
                (
                    FixedBitSet::with_capacity_and_blocks(
                        self.variables,
                        self.data.iter().map(|d| d[i]),
                    ),
                    *n,
                )
            })
    }

    fn to_icc(&self, icc: &FixedBitSet) -> Self {
        let mut new_icc_map = HashMap::with_hasher(S::default());
        for (o, n) in self.mask(icc).iter() {
            new_icc_map.entry(o).and_modify(|v| *v += n).or_insert(n);
        }
        Self::from(new_icc_map)
    }
}
