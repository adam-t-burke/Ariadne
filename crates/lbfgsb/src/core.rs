// SPDX-License-Identifier: BSD-3-Clause
//! Safe scalar line-search core translated from L-BFGS-B 3.0.
//!
//! Upstream: <http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>.
//! The upstream distribution's verbatim license and attribution are preserved
//! in `UPSTREAM_LICENSE.txt` and `THIRD_PARTY_NOTICES.md`; no upstream
//! copyright holder is inferred from its blank license template.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SearchWarning {
    RoundingErrors,
    StepTolerance,
    StepAtMaximum,
    StepAtMinimum,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SearchError {
    StepBelowMinimum,
    StepAboveMaximum,
    NonDescentDirection,
    NegativeFunctionTolerance,
    NegativeGradientTolerance,
    NegativeStepTolerance,
    NegativeMinimumStep,
    InvalidStepBounds,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum SearchResult {
    Evaluate(f64),
    Converged(f64),
    Warning { step: f64, reason: SearchWarning },
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct SearchConfig {
    pub(crate) function_tolerance: f64,
    pub(crate) gradient_tolerance: f64,
    pub(crate) step_tolerance: f64,
    pub(crate) minimum_step: f64,
    pub(crate) maximum_step: f64,
}

impl SearchConfig {
    pub(crate) fn lbfgsb(maximum_step: f64) -> Self {
        Self {
            function_tolerance: 1.0e-3,
            gradient_tolerance: 0.9,
            step_tolerance: 0.1,
            minimum_step: 0.0,
            maximum_step,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct StepInterval {
    pub(crate) best_step: f64,
    pub(crate) best_value: f64,
    pub(crate) best_derivative: f64,
    pub(crate) other_step: f64,
    pub(crate) other_value: f64,
    pub(crate) other_derivative: f64,
    pub(crate) bracketed: bool,
}

/// Exact four-case safeguarded step update from L-BFGS-B 3.0 `dcstep`.
pub(crate) fn dcstep(
    interval: &mut StepInterval,
    step: &mut f64,
    value: f64,
    derivative: f64,
    minimum_step: f64,
    maximum_step: f64,
) {
    const P66: f64 = 0.66;
    let stx = interval.best_step;
    let fx = interval.best_value;
    let dx = interval.best_derivative;
    let sty = interval.other_step;
    let fy = interval.other_value;
    let dy = interval.other_derivative;
    let stp = *step;
    let sgnd = derivative * (dx / dx.abs());

    let stpf = if value > fx {
        let theta = 3.0 * (fx - value) / (stp - stx) + dx + derivative;
        let s = theta.abs().max(dx.abs()).max(derivative.abs());
        let mut gamma = s * ((theta / s).powi(2) - (dx / s) * (derivative / s)).sqrt();
        if stp < stx {
            gamma = -gamma;
        }
        let p = (gamma - dx) + theta;
        let q = ((gamma - dx) + gamma) + derivative;
        let stpc = stx + (p / q) * (stp - stx);
        let stpq = stx + (dx / ((fx - value) / (stp - stx) + dx) / 2.0) * (stp - stx);
        interval.bracketed = true;
        if (stpc - stx).abs() < (stpq - stx).abs() {
            stpc
        } else {
            stpc + (stpq - stpc) / 2.0
        }
    } else if sgnd < 0.0 {
        let theta = 3.0 * (fx - value) / (stp - stx) + dx + derivative;
        let s = theta.abs().max(dx.abs()).max(derivative.abs());
        let mut gamma = s * ((theta / s).powi(2) - (dx / s) * (derivative / s)).sqrt();
        if stp > stx {
            gamma = -gamma;
        }
        let p = (gamma - derivative) + theta;
        let q = ((gamma - derivative) + gamma) + dx;
        let stpc = stp + (p / q) * (stx - stp);
        let stpq = stp + (derivative / (derivative - dx)) * (stx - stp);
        interval.bracketed = true;
        if (stpc - stp).abs() > (stpq - stp).abs() {
            stpc
        } else {
            stpq
        }
    } else if derivative.abs() < dx.abs() {
        let theta = 3.0 * (fx - value) / (stp - stx) + dx + derivative;
        let s = theta.abs().max(dx.abs()).max(derivative.abs());
        let mut gamma = s * 0.0_f64
            .max((theta / s).powi(2) - (dx / s) * (derivative / s))
            .sqrt();
        if stp > stx {
            gamma = -gamma;
        }
        let p = (gamma - derivative) + theta;
        let q = (gamma + (dx - derivative)) + gamma;
        let r = p / q;
        let stpc = if r < 0.0 && gamma != 0.0 {
            stp + r * (stx - stp)
        } else if stp > stx {
            maximum_step
        } else {
            minimum_step
        };
        let stpq = stp + (derivative / (derivative - dx)) * (stx - stp);
        if interval.bracketed {
            let candidate = if (stpc - stp).abs() < (stpq - stp).abs() {
                stpc
            } else {
                stpq
            };
            if stp > stx {
                candidate.min(stp + P66 * (sty - stp))
            } else {
                candidate.max(stp + P66 * (sty - stp))
            }
        } else {
            let candidate = if (stpc - stp).abs() > (stpq - stp).abs() {
                stpc
            } else {
                stpq
            };
            candidate.clamp(minimum_step, maximum_step)
        }
    } else if interval.bracketed {
        let theta = 3.0 * (value - fy) / (sty - stp) + dy + derivative;
        let s = theta.abs().max(dy.abs()).max(derivative.abs());
        let mut gamma = s * ((theta / s).powi(2) - (dy / s) * (derivative / s)).sqrt();
        if stp > sty {
            gamma = -gamma;
        }
        let p = (gamma - derivative) + theta;
        let q = ((gamma - derivative) + gamma) + dy;
        stp + (p / q) * (sty - stp)
    } else if stp > stx {
        maximum_step
    } else {
        minimum_step
    };

    if value > fx {
        interval.other_step = stp;
        interval.other_value = value;
        interval.other_derivative = derivative;
    } else {
        if sgnd < 0.0 {
            interval.other_step = stx;
            interval.other_value = fx;
            interval.other_derivative = dx;
        }
        interval.best_step = stp;
        interval.best_value = value;
        interval.best_derivative = derivative;
    }
    *step = stpf;
}

#[derive(Clone, Debug)]
pub(crate) struct DcSearch {
    config: SearchConfig,
    interval: StepInterval,
    stage: u8,
    initial_value: f64,
    initial_derivative: f64,
    test_derivative: f64,
    interval_minimum: f64,
    interval_maximum: f64,
    width: f64,
    previous_width: f64,
    current_step: f64,
}

impl DcSearch {
    pub(crate) fn start(
        value: f64,
        derivative: f64,
        step: f64,
        config: SearchConfig,
    ) -> Result<Self, SearchError> {
        if step < config.minimum_step {
            return Err(SearchError::StepBelowMinimum);
        }
        if step > config.maximum_step {
            return Err(SearchError::StepAboveMaximum);
        }
        if derivative >= 0.0 {
            return Err(SearchError::NonDescentDirection);
        }
        if config.function_tolerance < 0.0 {
            return Err(SearchError::NegativeFunctionTolerance);
        }
        if config.gradient_tolerance < 0.0 {
            return Err(SearchError::NegativeGradientTolerance);
        }
        if config.step_tolerance < 0.0 {
            return Err(SearchError::NegativeStepTolerance);
        }
        if config.minimum_step < 0.0 {
            return Err(SearchError::NegativeMinimumStep);
        }
        if config.maximum_step < config.minimum_step {
            return Err(SearchError::InvalidStepBounds);
        }
        let width = config.maximum_step - config.minimum_step;
        Ok(Self {
            config,
            interval: StepInterval {
                best_step: 0.0,
                best_value: value,
                best_derivative: derivative,
                other_step: 0.0,
                other_value: value,
                other_derivative: derivative,
                bracketed: false,
            },
            stage: 1,
            initial_value: value,
            initial_derivative: derivative,
            test_derivative: config.function_tolerance * derivative,
            interval_minimum: 0.0,
            interval_maximum: step + 4.0 * step,
            width,
            previous_width: width / 0.5,
            current_step: step,
        })
    }

    pub(crate) fn step(&self) -> f64 {
        self.current_step
    }

    pub(crate) fn advance(&mut self, value: f64, derivative: f64) -> SearchResult {
        let ftest = self.initial_value + self.current_step * self.test_derivative;
        if self.stage == 1 && value <= ftest && derivative >= 0.0 {
            self.stage = 2;
        }

        let mut warning = None;
        if self.interval.bracketed
            && (self.current_step <= self.interval_minimum
                || self.current_step >= self.interval_maximum)
        {
            warning = Some(SearchWarning::RoundingErrors);
        }
        if self.interval.bracketed
            && self.interval_maximum - self.interval_minimum
                <= self.config.step_tolerance * self.interval_maximum
        {
            warning = Some(SearchWarning::StepTolerance);
        }
        if self.current_step == self.config.maximum_step
            && value <= ftest
            && derivative <= self.test_derivative
        {
            warning = Some(SearchWarning::StepAtMaximum);
        }
        if self.current_step == self.config.minimum_step
            && (value > ftest || derivative >= self.test_derivative)
        {
            warning = Some(SearchWarning::StepAtMinimum);
        }

        if value <= ftest
            && derivative.abs() <= self.config.gradient_tolerance * (-self.initial_derivative)
        {
            return SearchResult::Converged(self.current_step);
        }
        if let Some(reason) = warning {
            return SearchResult::Warning {
                step: self.current_step,
                reason,
            };
        }

        if self.stage == 1 && value <= self.interval.best_value && value > ftest {
            let trial_step = self.current_step;
            let mut modified = StepInterval {
                best_step: self.interval.best_step,
                best_value: self.interval.best_value
                    - self.interval.best_step * self.test_derivative,
                best_derivative: self.interval.best_derivative - self.test_derivative,
                other_step: self.interval.other_step,
                other_value: self.interval.other_value
                    - self.interval.other_step * self.test_derivative,
                other_derivative: self.interval.other_derivative - self.test_derivative,
                bracketed: self.interval.bracketed,
            };
            dcstep(
                &mut modified,
                &mut self.current_step,
                value - trial_step * self.test_derivative,
                derivative - self.test_derivative,
                self.interval_minimum,
                self.interval_maximum,
            );
            self.interval = StepInterval {
                best_step: modified.best_step,
                best_value: modified.best_value + modified.best_step * self.test_derivative,
                best_derivative: modified.best_derivative + self.test_derivative,
                other_step: modified.other_step,
                other_value: modified.other_value + modified.other_step * self.test_derivative,
                other_derivative: modified.other_derivative + self.test_derivative,
                bracketed: modified.bracketed,
            };
        } else {
            dcstep(
                &mut self.interval,
                &mut self.current_step,
                value,
                derivative,
                self.interval_minimum,
                self.interval_maximum,
            );
        }

        if self.interval.bracketed {
            if (self.interval.other_step - self.interval.best_step).abs()
                >= 0.66 * self.previous_width
            {
                self.current_step = self.interval.best_step
                    + 0.5 * (self.interval.other_step - self.interval.best_step);
            }
            self.previous_width = self.width;
            self.width = (self.interval.other_step - self.interval.best_step).abs();
            self.interval_minimum = self.interval.best_step.min(self.interval.other_step);
            self.interval_maximum = self.interval.best_step.max(self.interval.other_step);
        } else {
            self.interval_minimum =
                self.current_step + 1.1 * (self.current_step - self.interval.best_step);
            self.interval_maximum =
                self.current_step + 4.0 * (self.current_step - self.interval.best_step);
        }
        self.current_step = self
            .current_step
            .clamp(self.config.minimum_step, self.config.maximum_step);
        if self.interval.bracketed
            && (self.current_step <= self.interval_minimum
                || self.current_step >= self.interval_maximum
                || self.interval_maximum - self.interval_minimum
                    <= self.config.step_tolerance * self.interval_maximum)
        {
            self.current_step = self.interval.best_step;
        }
        SearchResult::Evaluate(self.current_step)
    }
}

pub(crate) fn feasible_maximum_step(
    x: &[f64],
    direction: &[f64],
    lower: &[f64],
    upper: &[f64],
    iteration: usize,
) -> f64 {
    let constrained = lower
        .iter()
        .zip(upper)
        .any(|(&l, &u)| l.is_finite() || u.is_finite());
    if !constrained {
        return 1.0e10;
    }
    if iteration == 0 {
        return 1.0;
    }
    let mut maximum = 1.0e10;
    for (((&xi, &di), &lower), &upper) in x.iter().zip(direction).zip(lower).zip(upper) {
        if di < 0.0 && lower.is_finite() {
            let distance = lower - xi;
            if distance >= 0.0 {
                maximum = 0.0;
            } else if di * maximum < distance {
                maximum = distance / di;
            }
        } else if di > 0.0 && upper.is_finite() {
            let distance = upper - xi;
            if distance <= 0.0 {
                maximum = 0.0;
            } else if di * maximum > distance {
                maximum = distance / di;
            }
        }
    }
    maximum
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(bracketed: bool) -> StepInterval {
        StepInterval {
            best_step: 0.0,
            best_value: 0.0,
            best_derivative: -1.0,
            other_step: 2.0,
            other_value: 2.0,
            other_derivative: 1.0,
            bracketed,
        }
    }

    #[test]
    fn dcstep_exercises_all_four_fortran_cases() {
        let cases = [
            (false, 2.0, -0.5, true),    // higher value
            (false, -0.5, 0.5, true),    // opposite derivative
            (false, -0.5, -0.25, false), // decreased magnitude
            (false, -0.5, -2.0, false),  // non-decreased magnitude
        ];
        for (bracketed, value, derivative, expected_bracketed) in cases {
            let mut state = interval(bracketed);
            let mut step = 1.0;
            dcstep(&mut state, &mut step, value, derivative, 0.0, 4.0);
            assert!(step.is_finite());
            assert!((0.0..=4.0).contains(&step));
            assert_eq!(state.bracketed, expected_bracketed);
        }
        let mut state = interval(true);
        let mut step = 1.0;
        dcstep(&mut state, &mut step, -0.5, -2.0, 0.0, 4.0);
        assert!(step.is_finite(), "bracketed fourth case uses cubic step");
    }

    #[test]
    fn dcsrch_reports_convergence_and_boundary_warnings() {
        let mut search = DcSearch::start(0.0, -1.0, 1.0, SearchConfig::lbfgsb(1.0)).unwrap();
        assert_eq!(search.advance(-0.5, 0.0), SearchResult::Converged(1.0));

        let mut search = DcSearch::start(0.0, -1.0, 1.0, SearchConfig::lbfgsb(1.0)).unwrap();
        assert_eq!(
            search.advance(-2.0, -1.0),
            SearchResult::Warning {
                step: 1.0,
                reason: SearchWarning::StepAtMaximum
            }
        );

        let mut search = DcSearch::start(0.0, -1.0, 0.0, SearchConfig::lbfgsb(1.0)).unwrap();
        assert_eq!(
            search.advance(1.0, 0.0),
            SearchResult::Warning {
                step: 0.0,
                reason: SearchWarning::StepAtMinimum
            }
        );

        let mut search = DcSearch::start(0.0, -1.0, 1.0, SearchConfig::lbfgsb(2.0)).unwrap();
        search.interval.bracketed = true;
        search.interval_minimum = 1.0;
        search.interval_maximum = 1.5;
        assert_eq!(
            search.advance(1.0, -1.0),
            SearchResult::Warning {
                step: 1.0,
                reason: SearchWarning::RoundingErrors
            }
        );

        let mut search = DcSearch::start(0.0, -1.0, 1.0, SearchConfig::lbfgsb(2.0)).unwrap();
        search.interval.bracketed = true;
        search.interval_minimum = 0.96;
        search.interval_maximum = 1.04;
        assert_eq!(
            search.advance(1.0, -1.0),
            SearchResult::Warning {
                step: 1.0,
                reason: SearchWarning::StepTolerance
            }
        );
    }

    #[test]
    fn dcsrch_matches_official_reference_cases() {
        let fixture = include_str!("../tests/reference/fixtures/dcsrch-cases.csv");
        let records = fixture
            .lines()
            .skip(1)
            .map(|line| {
                let mut fields = line.splitn(3, ',');
                let name = fields.next().unwrap();
                let step = fields.next().unwrap().trim().parse::<f64>().unwrap();
                let task = fields.next().unwrap();
                (name, step, task)
            })
            .collect::<Vec<_>>();
        assert_eq!(
            records,
            [
                ("step-at-maximum", 1.0, "WARNING: STP = STPMAX"),
                ("step-at-minimum", 0.0, "WARNING: STP = STPMIN"),
                ("converged", 1.0, "CONVERGENCE"),
            ]
        );

        let mut maximum = DcSearch::start(0.0, -1.0, 1.0, SearchConfig::lbfgsb(1.0)).unwrap();
        assert_eq!(
            maximum.advance(-2.0, -1.0),
            SearchResult::Warning {
                step: records[0].1,
                reason: SearchWarning::StepAtMaximum,
            }
        );

        let mut minimum = DcSearch::start(0.0, -1.0, 0.0, SearchConfig::lbfgsb(1.0)).unwrap();
        assert_eq!(
            minimum.advance(1.0, 0.0),
            SearchResult::Warning {
                step: records[1].1,
                reason: SearchWarning::StepAtMinimum,
            }
        );

        let mut converged = DcSearch::start(0.0, -1.0, 1.0, SearchConfig::lbfgsb(1.0)).unwrap();
        assert_eq!(
            converged.advance(-0.5, 0.0),
            SearchResult::Converged(records[2].1)
        );
    }

    #[test]
    fn stpmx_handles_one_sided_and_unbounded_variables() {
        assert_eq!(
            feasible_maximum_step(
                &[1.0, 2.0],
                &[-2.0, 3.0],
                &[0.0, f64::NEG_INFINITY],
                &[f64::INFINITY, 5.0],
                1,
            ),
            0.5
        );
        assert_eq!(
            feasible_maximum_step(&[1.0], &[2.0], &[f64::NEG_INFINITY], &[f64::INFINITY], 4,),
            1.0e10
        );
    }
}
