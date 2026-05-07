#[derive(Debug, Clone, Copy)]
pub enum TempType {
    Exponential { per_step: f64 },
    Logarithmic,
    Linear { slope: f64 },
    Constant { steps: usize },
}

#[derive(Debug)]
pub struct AnnealingTemperature {
    pub(crate) start: f64,
    pub(crate) end: f64,
    pub(crate) temp_type: TempType,
    pub(crate) next_temp: Option<Box<AnnealingTemperature>>,
}

impl Default for AnnealingTemperature {
    fn default() -> Self {
        AnnealingTemperature {
            start: 500.0,
            end: 5.0,
            temp_type: TempType::Exponential { per_step: 0.1 },
            next_temp: None,
        }
    }
}

impl AnnealingTemperature {
    /// Create a new temperature curve for use in simulated annealing. The temperature
    /// will start at `start` and then be decayed exponentially by removing
    /// `temp * decrease_per_step` each step, until `temp < end`.
    ///
    /// The curve can be extended with an additional curve with different `end` and
    /// `decrease_per_step` values with the `then()` method.
    fn new(start: f64, end: f64, temp_type: TempType) -> AnnealingTemperature {
        AnnealingTemperature {
            start,
            end,
            temp_type,
            next_temp: None,
        }
    }

    /// Create a new exponentail temperature curve for use in simulated annealing.
    /// The temperature will start at `start` and then be decayed exponentially by
    /// subtracting `temp * per_step` each step, until `temp < end`.
    ///
    /// The curve can be extended with an additional curve with a different function
    /// and parameters with the `then_<function>()` methods.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature};
    /// let temp_curve = AnnealingTemperature::exponential(10_000.0, 5.0, 0.0003);
    /// ```
    pub fn exponential(start: f64, end: f64, per_step: f64) -> AnnealingTemperature {
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        assert!(
            per_step < 1.0 && per_step > 0.0,
            "per step should be between zero and one."
        );
        AnnealingTemperature::new(start, end, TempType::Exponential { per_step })
    }

    /// Create a new logarithmic curve for use in simulated annealing. The temperature
    /// will start at `start` and then decay logarithmically until `temp < end`.
    ///
    /// The curve can be extended with an additional curve with a different function
    /// and parameters with the `then_<function>()` methods.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature};
    /// let temp_curve = AnnealingTemperature::logarithmic(100_000.0, 0.05);
    /// ```
    pub fn logarithmic(start: f64, end: f64) -> AnnealingTemperature {
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        AnnealingTemperature::new(start, end, TempType::Logarithmic)
    }

    /// Create a new linear temperature curve for use in simulated annealing.
    /// The temperature will start at `start` and then decay linearly. by removing
    /// `slope` each step, until `temp < end`. `slope` should be larger than zero.
    ///
    /// The curve can be extended with an additional curve with a different function
    /// and parameters with the `then_<function>()` methods.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature};
    /// let temp_curve = AnnealingTemperature::linear(10_000.0, 5.0, 1.0);
    /// ```
    pub fn linear(start: f64, end: f64, slope: f64) -> AnnealingTemperature {
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        assert!(slope > 0.0, "Slope should be higher than zero.");
        AnnealingTemperature::new(start, end, TempType::Linear { slope })
    }

    /// Keeps the simulated annealing solver at a given temperature for the given
    /// amount of `steps`.
    ///
    /// The curve can be extended with an additional curve with a different function
    /// and parameters with the `then_<function>()` methods.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature};
    /// let temp_curve = AnnealingTemperature::constant(5.0, 10_000);
    /// ```
    pub fn constant(start: f64, steps: usize) -> AnnealingTemperature {
        assert!(steps > 0, "Steps should be greater than zero.");
        AnnealingTemperature::new(start, start, TempType::Constant { steps })
    }

    /// Counts and returns the amount of steps in the resulting iterator.
    pub fn steps(&self) -> usize {
        self.create_iter().count()
    }

    pub(crate) fn set_next_none(&mut self, next_value: AnnealingTemperature) {
        match &mut self.next_temp {
            Some(at) => at.set_next_none(next_value),
            None => self.next_temp = Some(Box::new(next_value)),
        };
    }

    pub(crate) fn get_last_end(&self) -> f64 {
        match &self.next_temp {
            Some(at) => at.get_last_end(),
            None => self.end,
        }
    }

    /// Extend this temperature curve with an additional exponential curve.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature, TempType};
    /// let temp_curve = AnnealingTemperature::exponential(10_000.0, 5.0, 0.0003)
    ///     .then_exponential(0.001, 0.00001);
    /// ```
    pub fn then_exponential(mut self, end: f64, per_step: f64) -> Self {
        let start = self.get_last_end();
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        assert!(
            per_step < 1.0 && per_step > 0.0,
            "per step should be between zero and one."
        );
        self.set_next_none(AnnealingTemperature::new(
            start,
            end,
            TempType::Exponential { per_step },
        ));
        self
    }

    /// Extend this temperature curve with an additional logarithmic curve.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature, TempType};
    /// let temp_curve = AnnealingTemperature::exponential(10_000.0, 5.0, 0.0003)
    ///     .then_logarithmic(0.1);
    /// ```
    pub fn then_logarithmic(mut self, end: f64) -> Self {
        let start = self.get_last_end();
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        self.set_next_none(AnnealingTemperature::new(start, end, TempType::Logarithmic));
        self
    }

    /// Extend this temperature curve with an additional linear curve.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature, TempType};
    /// let temp_curve = AnnealingTemperature::exponential(10_000.0, 5.0, 0.0003)
    ///     .then_linear(0.01, 0.1);
    /// ```
    pub fn then_linear(mut self, end: f64, slope: f64) -> Self {
        let start = self.get_last_end();
        assert!(end > 0.0, "End temperature should be greater than zero.");
        assert!(
            start > end,
            "Starting temperature should be greater than end temperature."
        );
        self.set_next_none(AnnealingTemperature::new(
            start,
            end,
            TempType::Linear { slope },
        ));
        self
    }

    /// Extend this temperature curve by keeping it at the current temperature for
    /// the given amount of `steps`.
    ///
    /// # Examples
    /// ```
    /// # use mcm_finder_lib::solvers::{AnnealingTemperature, TempType};
    /// let temp_curve = AnnealingTemperature::exponential(10_000.0, 5.0, 0.0003)
    ///     .then_constant(10_000);
    /// ```
    pub fn then_constant(mut self, steps: usize) -> Self {
        let start = self.get_last_end();
        self.set_next_none(AnnealingTemperature::new(
            start,
            start,
            TempType::Constant { steps },
        ));
        self
    }

    /// Returns iterator over the temperature values for this curve.
    pub fn create_iter(&self) -> AnnealTempIter {
        let mut iterator = AnnealTempIter::new(self.start, self.end, self.temp_type, None);
        iterator.next_temp = self.next_temp.as_ref().map(|t| Box::new(t.create_iter()));
        iterator
    }
}

#[derive(Debug, Clone)]
pub struct AnnealTempIter {
    start: f64,
    temp: f64,
    step: usize,
    end: f64,
    temp_type: TempType,
    next_temp: Option<Box<AnnealTempIter>>,
}

impl AnnealTempIter {
    pub(crate) fn new(
        temp: f64,
        end: f64,
        temp_type: TempType,
        next_temp: Option<Box<AnnealTempIter>>,
    ) -> AnnealTempIter {
        AnnealTempIter {
            start: temp,
            temp,
            step: 0,
            end,
            temp_type,
            next_temp,
        }
    }

    pub fn get_current_target(&self) -> &f64 {
        &self.end
    }

    fn delta(&mut self) {
        match self.temp_type {
            TempType::Exponential { per_step } => self.temp -= self.temp * per_step,
            TempType::Logarithmic => self.temp = self.start / (1.0 + self.step as f64),
            TempType::Linear { slope } => self.temp -= slope,
            TempType::Constant { steps } => {
                self.temp = self.start;
                if self.step == steps {
                    self.end = f64::INFINITY;
                };
            }
        }
    }
}

impl Iterator for AnnealTempIter {
    type Item = (f64, f64);

    fn next(&mut self) -> Option<Self::Item> {
        let output = if self.temp.total_cmp(&self.end).is_ge() {
            Some((self.temp, self.end))
        } else {
            match &self.next_temp {
                Some(t) => {
                    self.start = t.start;
                    self.end = t.end;
                    self.step = 0;
                    self.temp_type = t.temp_type;
                    self.next_temp = t.next_temp.clone();
                    Some((self.temp, self.end))
                }
                None => return None,
            }
        };
        self.delta();
        self.step += 1;
        output
    }
}
