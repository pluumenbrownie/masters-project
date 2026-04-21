mod mcm_error;

pub mod dataset;

pub mod mcm;

// #[derive(Debug, Default)]
// pub struct BasisSet {
//     basis_vectors: Vec<FixedBitSet>,
// }

// impl BasisSet {
//     pub fn new(vectors: Vec<FixedBitSet>) -> BasisSet {
//         BasisSet {
//             basis_vectors: vectors,
//         }
//     }

//     pub fn add(&mut self, basis: FixedBitSet) -> Result<(), MCMError> {
//         if self.basis_vectors.iter().any(|b| !b.is_disjoint(&basis)) {
//             return Err(MCMError::OverlappingBasis);
//         }
//         self.basis_vectors.push(basis);
//         Ok(())
//     }

// pub fn create_kset(&self, dataset: Dataset) -> KSet {
//     let new_vectors: Vec<(FixedBitSet, usize)> = dataset
//         .data
//         .iter()
//         .map(|(i, &n)| self.transform_mu_basis(i, n))
//         .collect();
//     // transformed mu vectors are not guaranteed to be unique, so we want to
//     // add them together without loosing information
//     let mut hashmap = HashMap::new();
//     for (o, n) in new_vectors {
//         hashmap.entry(o).and_modify(|v| *v += n).or_insert(n);
//     }
//     KSet {
//         data: hashmap,
//         datapoints: dataset.datapoints,
//     }
// }

// fn transform_mu_basis(&self, i: &FixedBitSet, n: usize) -> (FixedBitSet, usize) {
//     let mut o = FixedBitSet::with_capacity(self.basis_vectors.len());
//     for (nr, b) in self.basis_vectors.iter().enumerate() {
//         if (i & b).count_ones(..) % 2 == 1 {
//             o.insert(nr);
//         }
//     }
//     (o, n)
// }
// }

// #[derive(Debug)]
// pub struct KSet {
//     pub data: HashMap<FixedBitSet, usize>,
//     pub datapoints: usize,
// }
