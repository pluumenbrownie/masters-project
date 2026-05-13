use std::{cell::RefCell, collections::HashMap};

use fixedbitset::FixedBitSet;
use kdam::{Bar, BarExt, tqdm};
use rand::{
    RngExt,
    rngs::ThreadRng,
    seq::{IndexedRandom, IteratorRandom},
};

use crate::{
    dataset::{Dataset, VecDataset},
    mcm::MinimallyComplexModel,
    solvers::{Solver, SolverReport, get_log_e_cache},
};

/// The way individuals are chosen, proportional to their fitness.
pub enum SelectionType {
    /// Chance of being chosen is dependent linearly on individuals fitness.
    Linear,
    /// Chance of being chosen is dependent linearly on individuals squared fitness.
    Quadratic,
    /// Chance of being chosen is dependent exponentially on individuals fitness.
    Exponential,
    /// Take `from` individuals at random, select `choose` best ones.
    Tournament { choose: usize, from: usize },
}

impl SelectionType {
    fn select<'a>(
        &self,
        population: &[(&'a MinimallyComplexModel, f64)],
        amount: usize,
        rng: &mut ThreadRng,
    ) -> Vec<(&'a MinimallyComplexModel, f64)> {
        if self.is_ranked() {
            self.rank_selection(population, amount, rng)
        } else {
            self.tournament_selection(population, amount, rng)
        }
    }

    fn is_ranked(&self) -> bool {
        match &self {
            SelectionType::Linear => true,
            SelectionType::Quadratic => true,
            SelectionType::Exponential => true,
            SelectionType::Tournament { choose: _, from: _ } => false,
        }
    }

    fn rank_selection<'a>(
        &self,
        population: &[(&'a MinimallyComplexModel, f64)],
        amount: usize,
        rng: &mut ThreadRng,
    ) -> Vec<(&'a MinimallyComplexModel, f64)> {
        population
            .sample_weighted(rng, amount, |(_mcm, w)| match &self {
                SelectionType::Linear => *w,
                SelectionType::Quadratic => w.powi(2),
                SelectionType::Exponential => 1.0 - (-w).exp(),
                SelectionType::Tournament { choose: _, from: _ } => {
                    unreachable!("Tournament selection should have been eliminated earlier.")
                }
            })
            .unwrap()
            .cloned()
            .collect()
    }

    /// Use tournament selection to select `amount` unique individuals.
    ///
    /// This function is really gross to ensure the chosen individuals are
    /// all unique.
    fn tournament_selection<'a>(
        &self,
        population: &[(&'a MinimallyComplexModel, f64)],
        amount: usize,
        rng: &mut ThreadRng,
    ) -> Vec<(&'a MinimallyComplexModel, f64)> {
        match &self {
            SelectionType::Tournament { choose, from } => {
                let mut output = vec![];
                let mut picks = 0usize;
                let mut population = population.to_vec();

                while picks > 0 {
                    let amount_round = *choose.max(&(amount - picks));

                    let tournament_indexes = (0..population.len()).sample(rng, *from);
                    let mut tournament: Vec<_> = tournament_indexes
                        .into_iter()
                        .map(|i| {
                            let pair = population[i];
                            (pair.0, pair.1, i)
                        })
                        .collect();
                    tournament.sort_by(|a, b| b.1.total_cmp(&a.1));
                    let winners: Vec<_> = tournament.into_iter().take(amount_round).collect();

                    let mut winner_indexes: Vec<_> = winners.iter().map(|t| t.2).collect();
                    winner_indexes.sort();
                    for x in winner_indexes.into_iter().rev() {
                        population.remove(x);
                    }

                    output.extend(winners.into_iter().map(|(a, b, _)| (a, b)));

                    picks += amount_round;
                }

                output
            }
            _ => unreachable!("Ranked selection should have been eliminated earlier."),
        }
    }
}

pub struct EvolutionarySolver {
    dataset: VecDataset,
    generation_size: usize,
    generations: usize,
    shuffle_steps: usize,
    crossover_probability: f64,
    mutation_rate: usize,
    survivors: usize,
    parent_selection_type: SelectionType,
    survivor_selection_type: SelectionType,
    elitism: usize,
    rng: RefCell<ThreadRng>,
}

impl Solver for EvolutionarySolver {
    fn from_file(filepath: &std::path::Path) -> Result<Self, crate::mcm_error::MCMError>
    where
        Self: Sized,
    {
        let dataset = VecDataset::read_from_file(filepath)?;
        let variables = dataset.variables();
        Ok(EvolutionarySolver {
            dataset,
            generations: 10_000,
            generation_size: 4,
            shuffle_steps: 2 * variables,
            crossover_probability: 0.3,
            mutation_rate: 1,
            survivors: 2,
            parent_selection_type: SelectionType::Linear,
            survivor_selection_type: SelectionType::Linear,
            elitism: 0,
            rng: RefCell::new(rand::rng()),
        })
    }

    fn solve(&self) -> SolverReport {
        let mut log_e_cache = get_log_e_cache();
        let mut progress = tqdm!();

        let mut generation = self.generate_starting_population();

        for gen_nr in 0..self.generations {
            let fitness_generation =
                self.calculate_fitness(&mut log_e_cache, &mut progress, &generation);
            progress.set_description(format!(
                "Generation {} - Best Log E: {:.0}",
                gen_nr,
                fitness_generation
                    .iter()
                    .map(|x| x.0.log_e(&self.dataset, &mut log_e_cache))
                    .max_by(|a, b| a.total_cmp(b))
                    .unwrap()
            ));

            let survivors = self.survivor_selection(fitness_generation);
            generation = self.generate_next_generation(survivors);
        }

        let best_mcm = generation
            .iter()
            .map(|mcm| (mcm, mcm.log_e(&self.dataset, &mut log_e_cache)))
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        SolverReport::new(
            best_mcm.0.clone(),
            best_mcm.1,
            HashMap::from([(
                "Unique ICCs covered".into(),
                format!("{}", log_e_cache.unwrap().len()),
            )]),
        )
    }
}

impl EvolutionarySolver {
    fn generate_starting_population(&self) -> Vec<MinimallyComplexModel> {
        let mut output = Vec::with_capacity(self.generation_size);

        for _ in 0..self.generation_size {
            let mut mcm = MinimallyComplexModel::full(self.dataset.variables().try_into().unwrap());
            for _ in 0..self.shuffle_steps {
                mcm = mcm.mutate(&mut self.rng.borrow_mut());
            }
            output.push(mcm);
        }

        output
    }

    fn survivor_selection<'a>(
        &self,
        generation: Vec<(&'a MinimallyComplexModel, f64)>,
    ) -> Vec<(&'a MinimallyComplexModel, f64)> {
        let mut generation = generation;
        // elitism
        generation.sort_by(|a, b| a.1.total_cmp(&b.1));
        let mut output = generation.split_off(generation.len() - self.elitism);
        // normal selection
        output.extend(self.survivor_selection_type.select(
            &generation,
            self.survivors,
            &mut self.rng.borrow_mut(),
        ));
        output
    }

    fn generate_next_generation(
        &self,
        survivors: Vec<(&MinimallyComplexModel, f64)>,
    ) -> Vec<MinimallyComplexModel> {
        self.generate_children(survivors)
    }

    fn generate_children(
        &self,
        survivors: Vec<(&MinimallyComplexModel, f64)>,
    ) -> Vec<MinimallyComplexModel> {
        let mut next_generation: Vec<MinimallyComplexModel> = survivors
            .iter()
            .take(self.elitism)
            .map(|(mcm, _f)| *mcm)
            .cloned()
            .collect();

        let parent_couples = (self.generation_size - self.elitism).div_ceil(2);
        let parents = self.parent_selection(parent_couples, survivors);

        for (parent_a, parent_b) in parents {
            let mut child_a = parent_a.clone();
            let mut child_b = parent_b.clone();

            for _ in 0..self.mutation_rate {
                child_a = child_a.mutate(&mut self.rng.borrow_mut());
                child_b = child_b.mutate(&mut self.rng.borrow_mut());
            }

            self.crossover(&mut child_a, &mut child_b);

            next_generation.push(child_a);
            next_generation.push(child_b);
        }

        // next_generation could be too long if self.generation_size - self.elitism is uneven.
        if next_generation.len() > self.generation_size {
            next_generation.pop();
        }

        next_generation
    }

    fn parent_selection<'a>(
        &self,
        parent_couples: usize,
        survivors: Vec<(&'a MinimallyComplexModel, f64)>,
    ) -> Vec<(&'a MinimallyComplexModel, &'a MinimallyComplexModel)> {
        let mut output = vec![];

        for _ in 0..parent_couples {
            let pair = self
                .parent_selection_type
                .select(&survivors, 2, &mut self.rng.borrow_mut());
            output.push((pair[0].0, pair[1].0));
        }

        output
    }

    fn crossover(&self, a: &mut MinimallyComplexModel, b: &mut MinimallyComplexModel) {
        let zipped_up_vectors = a.to_vector().into_iter().zip(b.to_vector());
        let (a_vector, b_vector) = zipped_up_vectors
            .map(|(a, b)| {
                if self
                    .rng
                    .borrow_mut()
                    .random_bool(self.crossover_probability)
                {
                    (b, a)
                } else {
                    (a, b)
                }
            })
            .unzip();

        *a = MinimallyComplexModel::from_vector(a_vector);
        *b = MinimallyComplexModel::from_vector(b_vector);
    }

    fn calculate_fitness<'a>(
        &self,
        log_e_cache: &mut Option<HashMap<FixedBitSet, f64>>,
        progress: &mut Bar,
        generation: &'a [MinimallyComplexModel],
    ) -> Vec<(&'a MinimallyComplexModel, f64)> {
        let log_e: Vec<_> = generation
            .iter()
            .map(|mcm| {
                let _ = progress.update(1);
                (mcm, mcm.log_e(&self.dataset, log_e_cache))
            })
            .collect();
        let lowest = log_e.iter().min_by(|a, b| a.1.total_cmp(&b.1)).unwrap().1;
        let highest = log_e.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap().1;
        let delta = highest - lowest;
        log_e
            .into_iter()
            .map(|(mcm, log_e)| (mcm, (log_e - lowest + 1.0) / delta.abs()))
            .collect()
    }

    /// Set the maximum amount of generations run by this solver.
    pub fn set_generations(mut self, generations: usize) -> Self {
        self.generations = generations;
        self
    }

    /// Set the amount of indivituals in the population per generation.
    pub fn set_generation_size(mut self, size: usize) -> Self {
        self.generation_size = size;
        self
    }

    /// Set the amount of indivituals which survive after each generation.
    pub fn set_survivors(mut self, survivors: usize) -> Self {
        self.survivors = survivors;
        self
    }

    /// Set how many times the newly generated, random MCMs will be shuffled at the
    /// start of the algorithm.
    pub fn set_shuffle_steps(mut self, steps: usize) -> Self {
        self.shuffle_steps = steps;
        self
    }

    /// Set the amount of times the mutation function will be called on each new child.
    pub fn set_mutation_rate(mut self, rate: usize) -> Self {
        self.mutation_rate = rate;
        self
    }

    /// Set the crossover chance of the variables.
    pub fn set_crossover_probability(mut self, probability: f64) -> Self {
        self.crossover_probability = probability;
        self
    }

    /// Set the type of selection to use for parent selection.
    pub fn set_parent_selection(mut self, selection_type: SelectionType) -> Self {
        self.parent_selection_type = selection_type;
        self
    }

    /// Set the type of selection to use for survivor selection.
    pub fn set_survivor_selection(mut self, selection_type: SelectionType) -> Self {
        self.survivor_selection_type = selection_type;
        self
    }

    /// Set how many of the best generated models each generation are guaranteed to
    /// survive to the next generation.
    pub fn set_elitism(mut self, amount: usize) -> Self {
        self.elitism = amount;
        self
    }
}
