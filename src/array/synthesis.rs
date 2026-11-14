//! Array pattern synthesis: element excitations that keep the gain between
//! lower and upper bounds at sampled directions (contoured / shaped beams,
//! side-lobe regions, not-to-exceed regions), and closed-form nulls.
//!
//! # Method
//! [`synthesize`] uses alternating projections (Bucci, Franceschetti,
//! Mazzarella, Panariello, "Intersection approach to array pattern
//! synthesis", IEE Proc. H 137(6), 1990). The far field of the array at
//! sample `m` is `F_m = e_m Σ_n w_n exp(j(k₀ r_n·û_m + φ_steer,n))`, with
//! `e_m` the element field amplitude, and its gain `G_m = |F_m|² / Σ|w_n|²`
//! (the normalisation of [`PhasedArray`], ADR 0004). Each iteration:
//! 1. projects the field onto the bounds: `|F_m|` is clamped to
//!    `[√(G_min P), √(G_max P)]`, phase kept;
//! 2. moves the excitation towards that field: a few conjugate-gradient
//!    steps (CGLS) on `min ‖A w − F'‖²` from the current excitation;
//! 3. projects the excitation onto the hardware: failed elements off,
//!    [`ExcitationLimits`], phase quantisation of the total phase.
//!
//! The result is a local optimum that depends on the start (the array's
//! current excitation: warm start); the best iterate is returned. Bounds
//! may be relative to a *coverage level* that is either fixed or maximised
//! ([`CoverageLevel::Maximise`]): the level is then refitted at every
//! iteration while the shape converges (half the budget: aimed at the
//! median coverage gain, so dips are lifted), then raised in fixed-level runs
//! (1 dB steps, then bisection) as long as the bounds can still be met; if
//! the shaping phase does not meet them, fixed levels at and a few steps
//! below its own are tried first.
//!
//! Design is done at the design frequency, without the array's random
//! errors. Cost per iteration: about `10 · samples · elements` complex
//! exponentials (rows in parallel with the `parallel` feature).

use crate::array::PhasedArray;
use crate::array::excitation::{ExcitationLimits, quantise_phase};
use crate::direction::Direction;
use crate::error::{AntennaError, Result, invalid};
use crate::math::Complex;
use crate::pattern::GainPattern;

/// A gain level: absolute, or relative to the coverage level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Level {
    /// Absolute gain (dBi).
    Dbi(f64),
    /// Coverage level plus this offset (dB; negative = below).
    Coverage(f64),
}

/// Gain bounds at one direction (antenna frame).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PatternBound {
    /// Direction of the sample.
    pub direction: Direction,
    /// Lowest gain allowed (`None` = no lower bound).
    pub min: Option<Level>,
    /// Highest gain allowed (`None` = no upper bound).
    pub max: Option<Level>,
}

/// How the coverage level of [`Level::Coverage`] bounds is set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoverageLevel {
    /// A fixed level (dBi).
    Fixed(f64),
    /// The highest level that can be met within the tolerance (adaptive).
    Maximise,
}

/// Synthesis options.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SynthesisOptions {
    /// Coverage level of the relative bounds.
    pub coverage: CoverageLevel,
    /// Iteration budget (≥ 1).
    pub iterations: usize,
    /// Largest violation (dB) for a design to count as meeting the bounds.
    pub tolerance_db: f64,
    /// Amplitude limits of the beamformer.
    pub limits: ExcitationLimits,
}

impl Default for SynthesisOptions {
    fn default() -> Self {
        Self {
            coverage: CoverageLevel::Maximise,
            iterations: 300,
            tolerance_db: 0.5,
            limits: ExcitationLimits::default(),
        }
    }
}

/// Result of [`synthesize`].
#[derive(Debug, Clone, PartialEq)]
pub struct Synthesis {
    /// Element amplitudes, relative to the largest (failed elements 0).
    pub amplitudes: Vec<f64>,
    /// Element phase offsets (rad), added to the steering phase.
    pub phases_rad: Vec<f64>,
    /// Coverage level used (dBi); `None` without coverage-relative bounds.
    pub coverage_level_dbi: Option<f64>,
    /// Gain at each bound direction (dBi), in input order.
    pub gains_dbi: Vec<f64>,
    /// Largest violation of a bound (dB, 0 when all are met).
    pub max_violation_db: f64,
    /// `max_violation_db ≤ tolerance_db`.
    pub feasible: bool,
    /// Largest violation after each iteration (dB).
    pub history_db: Vec<f64>,
}

/// Fixed levels tried below the shaped one when it does not meet the bounds (dB).
const DOWN_STEPS_DB: [f64; 5] = [0.0, 1.0, 2.0, 4.0, 8.0];
/// First step when raising the coverage level (dB).
const LEVEL_STEP_DB: f64 = 1.0;
/// The coverage level is raised until it is known within this (dB).
const LEVEL_RESOLUTION_DB: f64 = 0.1;
/// Conjugate-gradient steps of the least-squares back-projection per iteration.
const CG_STEPS: usize = 4;
/// Smallest linear gain used for violations of a zero field.
const MIN_GAIN_LIN: f64 = 1e-30;

/// Excitations of `array` that keep its gain within `bounds`.
///
/// Starts from the array's current excitation (amplitudes and phase
/// offsets), keeps its failed elements off and its phase quantisation.
///
/// ```
/// use antenna_r::{Direction, GainPattern, Isotropic};
/// use antenna_r::array::{PhasedArray, lattice::linear};
/// use antenna_r::array::synthesis::{synthesize, CoverageLevel, Level, PatternBound, SynthesisOptions};
/// let lambda = rf_r::fspl::wavelength_m(10e9);
/// let a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &linear(16, lambda / 2.0).unwrap(), 10e9).unwrap();
/// // Keep the region beyond 15° at least 25 dB below the beam on boresight.
/// let mut bounds = vec![PatternBound { direction: Direction::BORESIGHT, min: Some(Level::Coverage(0.0)), max: None }];
/// for i in 0..60 {
///     let t = (15.0 + i as f64 * 1.25).to_radians();
///     for phi in [0.0, std::f64::consts::PI] {
///         bounds.push(PatternBound { direction: Direction::from_theta_phi(t, phi), min: None, max: Some(Level::Coverage(-25.0)) });
///     }
/// }
/// let opts = SynthesisOptions { coverage: CoverageLevel::Fixed(11.0), iterations: 400, ..Default::default() };
/// let s = synthesize(&a, &bounds, &opts).unwrap();
/// assert!(s.feasible, "violation {} dB", s.max_violation_db);
/// ```
///
/// # Errors
/// [`AntennaError::InvalidParameter`] for no bounds, a non-finite direction
/// or level, `iterations == 0`, a non-positive tolerance, or invalid limits.
pub fn synthesize<E: GainPattern>(
    array: &PhasedArray<E>,
    bounds: &[PatternBound],
    options: &SynthesisOptions,
) -> Result<Synthesis> {
    check(bounds, options)?;
    let problem = Problem::new(array, bounds);
    let (amps, phases) = array.excitation();
    let start: Vec<Complex> = amps
        .iter()
        .zip(phases)
        .map(|(a, p)| Complex::from_polar(*a, *p))
        .collect();
    Ok(problem.solve(bounds, &start, options))
}

/// Closed-form nulls: the excitation closest to the array's current one (in
/// the least-squares sense) whose array factor is zero towards every
/// direction of `nulls`:
/// `w = w₀ − Cᴴ(CCᴴ)⁻¹C w₀` with `C_kn = exp(j(k₀ r_n·û_k + φ_steer,n))`
/// (linearly constrained minimum change). Needs amplitude control
/// (amplitudes change); failed elements stay off.
///
/// ```
/// use antenna_r::{Direction, GainPattern, Isotropic};
/// use antenna_r::array::{PhasedArray, lattice::linear};
/// use antenna_r::array::synthesis::place_nulls;
/// let lambda = rf_r::fspl::wavelength_m(10e9);
/// let a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &linear(16, lambda / 2.0).unwrap(), 10e9).unwrap();
/// let null = Direction::from_theta_phi(20f64.to_radians(), 0.0);
/// let (amps, phases) = place_nulls(&a, &[null]).unwrap();
/// let b = a.with_excitation(&amps, &phases).unwrap().with_floor_dbi(-300.0).unwrap();
/// assert!(b.gain_dbi(null) < -100.0);
/// assert!(b.gain_dbi(Direction::BORESIGHT) > 11.5);
/// ```
///
/// # Errors
/// [`AntennaError::InvalidParameter`] for no null, a non-finite direction,
/// more nulls than radiating elements, nulls that cannot be met together
/// (e.g. duplicates) or that would remove the whole excitation.
pub fn place_nulls<E: GainPattern>(array: &PhasedArray<E>, nulls: &[Direction]) -> Result<(Vec<f64>, Vec<f64>)> {
    let dirs: Vec<[f64; 3]> = nulls.iter().map(|d| [d.u(), d.v(), d.w()]).collect();
    let (amps, phases) = array.excitation();
    let w0: Vec<Complex> = amps
        .iter()
        .zip(phases)
        .map(|(a, p)| Complex::from_polar(*a, *p))
        .collect();
    let geo = Geometry::of(array);
    null_steer(&geo, &dirs, &w0)
}

fn check(bounds: &[PatternBound], options: &SynthesisOptions) -> Result<()> {
    if bounds.is_empty() {
        return Err(invalid("bounds", "give at least one bound"));
    }
    let finite = |l: Option<Level>| match l {
        None => true,
        Some(Level::Dbi(v) | Level::Coverage(v)) => v.is_finite(),
    };
    for b in bounds {
        if b.direction.w().is_nan() {
            return Err(invalid("bounds", "directions must be finite"));
        }
        if !finite(b.min) || !finite(b.max) {
            return Err(invalid("bounds", "levels must be finite"));
        }
    }
    if options.iterations == 0 {
        return Err(invalid("iterations", "must be ≥ 1"));
    }
    if !(options.tolerance_db.is_finite() && options.tolerance_db > 0.0) {
        return Err(invalid("tolerance_db", "must be finite and > 0"));
    }
    if let CoverageLevel::Fixed(l) = options.coverage
        && !l.is_finite()
    {
        return Err(invalid("coverage", "the level must be finite"));
    }
    options.limits.validate()
}

/// Element geometry and steering of an array (non-generic).
struct Geometry {
    /// `k₀ x`, `k₀ y`, `k₀ z` per element.
    kpos: Vec<[f64; 3]>,
    /// Steering phase per element (rad).
    steer: Vec<f64>,
    /// Elements that radiate.
    on: Vec<bool>,
    /// Phase quantisation of the array.
    phase_bits: Option<u8>,
}

impl Geometry {
    fn of<E: GainPattern>(array: &PhasedArray<E>) -> Self {
        let k0 = std::f64::consts::TAU * array.design_freq_hz() / rf_r::constants::SPEED_OF_LIGHT;
        Self {
            kpos: array
                .positions()
                .iter()
                .map(|p| [k0 * p.x, k0 * p.y, k0 * p.z])
                .collect(),
            steer: array.steering_phases_rad(),
            on: array.failed().iter().map(|f| !f).collect(),
            phase_bits: array.phase_bits(),
        }
    }

    /// `exp(j(k₀ r_n·û + φ_steer,n))`.
    #[inline]
    fn term(&self, n: usize, d: &[f64; 3]) -> Complex {
        let p = &self.kpos[n];
        let phase = p[0].mul_add(d[0], p[1].mul_add(d[1], p[2] * d[2])) + self.steer[n];
        let (s, c) = phase.sin_cos();
        Complex::new(c, s)
    }
}

/// A synthesis problem: geometry plus samples.
struct Problem {
    geo: Geometry,
    /// Sample directions `(u, v, w)`.
    dirs: Vec<[f64; 3]>,
    /// Element field amplitude at each sample.
    elem: Vec<f64>,
}

impl Problem {
    fn new<E: GainPattern>(array: &PhasedArray<E>, bounds: &[PatternBound]) -> Self {
        let element = array.element();
        Self {
            geo: Geometry::of(array),
            dirs: bounds
                .iter()
                .map(|b| [b.direction.u(), b.direction.v(), b.direction.w()])
                .collect(),
            elem: bounds
                .iter()
                .map(|b| 10f64.powf(element.gain_dbi(b.direction) / 20.0))
                .collect(),
        }
    }

    /// `F = A w` (field at every sample).
    fn field(&self, w: &[Complex]) -> Vec<Complex> {
        rows(self.dirs.len(), |m| {
            let d = &self.dirs[m];
            let mut f = Complex::ZERO;
            for (n, wn) in w.iter().enumerate() {
                if self.geo.on[n] {
                    f += *wn * self.geo.term(n, d);
                }
            }
            f.scale(self.elem[m])
        })
    }

    /// `Aᴴ r` (back-projection onto the elements).
    fn adjoint(&self, r: &[Complex]) -> Vec<Complex> {
        rows(self.geo.kpos.len(), |n| {
            if !self.geo.on[n] {
                return Complex::ZERO;
            }
            let mut g = Complex::ZERO;
            for (m, rm) in r.iter().enumerate() {
                g += *rm * self.geo.term(n, &self.dirs[m]).conj().scale(self.elem[m]);
            }
            g
        })
    }

    fn solve(&self, bounds: &[PatternBound], start: &[Complex], options: &SynthesisOptions) -> Synthesis {
        let relative = bounds
            .iter()
            .any(|b| matches!(b.min, Some(Level::Coverage(_))) || matches!(b.max, Some(Level::Coverage(_))));
        let mode = match (options.coverage, relative) {
            (_, false) => Mode::Fixed(None),
            (CoverageLevel::Fixed(l), true) => Mode::Fixed(Some(l)),
            (CoverageLevel::Maximise, true) if has_coverage_min(bounds) => Mode::Adaptive,
            // Only relative upper bounds: the start's maximum gain is the level.
            (CoverageLevel::Maximise, true) => Mode::Fixed(Some(self.max_gain(start))),
        };
        let ctx = Run { bounds, options };
        let mut history = Vec::with_capacity(options.iterations);
        let start = self.project(start, options);
        let best = match mode {
            Mode::Adaptive => {
                // Shape at the natural level, then raise the level while it can be met.
                let shaped = (options.iterations / 2).max(1);
                let first = self.iterate(&ctx, start, mode, shaped, &mut history);
                self.raise(&ctx, first, options.iterations - shaped, &mut history)
            }
            Mode::Fixed(_) => self.iterate(&ctx, start, mode, options.iterations, &mut history),
        };
        let max_amp = best.w.iter().map(|c| c.abs()).fold(0.0, f64::max);
        Synthesis {
            amplitudes: best
                .w
                .iter()
                .map(|c| if max_amp > 0.0 { c.abs() / max_amp } else { 0.0 })
                .collect(),
            phases_rad: best.w.iter().map(|c| c.arg()).collect(),
            coverage_level_dbi: best.level,
            gains_dbi: best.gains,
            max_violation_db: best.violation,
            feasible: best.violation <= options.tolerance_db,
            history_db: history,
        }
    }

    /// Raises the coverage level after the shaping phase: fixed-level runs
    /// (1 dB steps, then bisection to 0.1 dB), each warm-started from the
    /// last feasible design and judged at the level it reaches (its weakest
    /// coverage sample). A shaping result that does not meet the bounds is
    /// first retried at fixed levels at and up to 8 dB below its own; if none
    /// meets them, it is shaped further with the remaining budget.
    fn raise(&self, ctx: &Run<'_>, first: Iterate, budget: usize, history: &mut Vec<f64>) -> Iterate {
        let (bounds, options) = (ctx.bounds, ctx.options);
        let chunk = (budget / 10).max(1);
        let mut used = 0;
        // Judge a fixed-level run at the level it reached.
        let run = |w: Vec<Complex>, trial: f64, n: usize, history: &mut Vec<f64>| {
            let mut r = self.iterate(ctx, w, Mode::Fixed(Some(trial)), n, history);
            let reached = level_bracket(bounds, &r.gains).1;
            r.violation = violation_at(bounds, &r.gains, Some(reached));
            r.level = Some(reached);
            r
        };
        let mut lo = first;
        if lo.violation > options.tolerance_db {
            // Fixed-level runs at and below the shaped level, a few steps
            // only; failing that, more shaping.
            let base = lo.level.unwrap_or(0.0);
            // Each step takes one chunk; at most `budget / chunk` of them fit.
            let steps = budget / chunk;
            let found = DOWN_STEPS_DB.iter().take(steps).find_map(|down| {
                used += chunk;
                let r = run(lo.w.clone(), base - down, chunk, history);
                (r.violation <= options.tolerance_db).then_some(r)
            });
            let rest = budget - used;
            lo = found.unwrap_or_else(|| {
                better(
                    self.iterate(ctx, lo.w.clone(), Mode::Adaptive, rest, history),
                    lo,
                    options.tolerance_db,
                )
            });
        }
        if lo.violation > options.tolerance_db {
            return lo;
        }
        let mut hi: Option<f64> = None;
        while used < budget {
            let l = lo.level.unwrap_or(0.0);
            let trial = hi.map_or(l + LEVEL_STEP_DB, |h| f64::midpoint(l, h));
            let n = chunk.min(budget - used);
            used += n;
            let r = run(lo.w.clone(), trial, n, history);
            if r.violation <= options.tolerance_db && r.level.unwrap_or(f64::NEG_INFINITY) > l {
                lo = r;
            } else {
                hi = Some(trial);
            }
            if hi.is_some_and(|h| h - lo.level.unwrap_or(0.0) < LEVEL_RESOLUTION_DB) {
                break;
            }
        }
        lo
    }

    /// Highest gain (dBi) over the samples.
    fn max_gain(&self, w: &[Complex]) -> f64 {
        let p = power(w, &self.geo.on);
        let peak = self.field(w).iter().map(|f| f.norm_sqr()).fold(0.0, f64::max);
        10.0 * (peak / p).max(MIN_GAIN_LIN).log10()
    }

    /// Projection iterations; the best iterate (meets the bounds at the
    /// highest level, else has the smallest violation).
    fn iterate(
        &self,
        ctx: &Run<'_>,
        w: Vec<Complex>,
        mode: Mode,
        iterations: usize,
        history: &mut Vec<f64>,
    ) -> Iterate {
        let mut current = self.evaluate(ctx.bounds, w, mode);
        let mut best = current.it.clone();
        for _ in 0..iterations {
            history.push(current.it.violation);
            if current.it.violation == 0.0 {
                break;
            }
            let next = self.step(ctx, &current);
            current = self.evaluate(ctx.bounds, next, mode);
            if current.it.better_than(&best, ctx.options.tolerance_db) {
                best = current.it.clone();
            }
        }
        best
    }

    /// Field, gains, levels and violation of an excitation.
    fn evaluate(&self, bounds: &[PatternBound], w: Vec<Complex>, mode: Mode) -> Evaluated {
        let field = self.field(&w);
        let power = power(&w, &self.geo.on);
        let gains: Vec<f64> = field
            .iter()
            .map(|fm| 10.0 * (fm.norm_sqr() / power).max(MIN_GAIN_LIN).log10())
            .collect();
        // Level the design is judged at (the weakest coverage sample), and
        // level it aims for: the median coverage sample, so dips are lifted,
        // or higher when the ripple is exceeded.
        let (level, target) = match mode {
            Mode::Fixed(l) => (l, l),
            Mode::Adaptive => {
                let (lo, hi) = level_bracket(bounds, &gains);
                (Some(hi), Some(median_coverage(bounds, &gains).max(lo)))
            }
        };
        let violation = violation_at(bounds, &gains, level);
        Evaluated {
            it: Iterate {
                w,
                level,
                violation,
                gains,
            },
            field,
            power,
            target,
        }
    }

    /// One projection step: field clamped to the bounds, gradient step,
    /// hardware projection.
    fn step(&self, ctx: &Run<'_>, e: &Evaluated) -> Vec<Complex> {
        let (bounds, options) = (ctx.bounds, ctx.options);
        // Lower bounds aim at the target, upper bounds at the judged level.
        let (c_lo, c_hi) = (e.target.unwrap_or(0.0), e.it.level.unwrap_or(0.0));
        let p = e.power;
        let r: Vec<Complex> = e
            .field
            .iter()
            .zip(bounds)
            .map(|(fm, b)| {
                let g = (fm.norm_sqr() / p).max(MIN_GAIN_LIN);
                let lo = linear(b.min, c_lo, 0.0);
                let hi = linear(b.max, c_hi, f64::INFINITY);
                let t = if g < lo {
                    (lo * p).sqrt()
                } else if g > hi {
                    (hi * p).sqrt()
                } else {
                    return Complex::ZERO;
                };
                // Same phase, magnitude t (an exactly zero field gets no correction).
                fm.scale(t / fm.abs().max(MIN_GAIN_LIN)) - *fm
            })
            .collect();
        let next = self.least_squares(&e.it.w, r);
        self.project(&next, options)
    }

    /// `CG_STEPS` conjugate-gradient steps on `min ‖A w − (A w₀ + r)‖²`
    /// from `w₀` (CGLS; `r` is the residual of the projected field).
    fn least_squares(&self, w0: &[Complex], mut r: Vec<Complex>) -> Vec<Complex> {
        let mut w = w0.to_vec();
        let mut s = self.adjoint(&r);
        let mut p = s.clone();
        let mut gamma = norm_sqr(&s);
        for _ in 0..CG_STEPS {
            let q = self.field(&p);
            let qq = norm_sqr(&q);
            if gamma <= MIN_GAIN_LIN || qq <= MIN_GAIN_LIN {
                break;
            }
            let alpha = gamma / qq;
            for (wn, pn) in w.iter_mut().zip(&p) {
                *wn += pn.scale(alpha);
            }
            for (rm, qm) in r.iter_mut().zip(&q) {
                *rm -= qm.scale(alpha);
            }
            s = self.adjoint(&r);
            let next = norm_sqr(&s);
            let beta = next / gamma;
            for (pn, sn) in p.iter_mut().zip(&s) {
                *pn = *sn + pn.scale(beta);
            }
            gamma = next;
        }
        w
    }

    /// Hardware projection: failed elements off, amplitude limits, total-phase quantisation.
    fn project(&self, w: &[Complex], options: &SynthesisOptions) -> Vec<Complex> {
        let amps: Vec<f64> = w
            .iter()
            .zip(&self.geo.on)
            .map(|(c, &on)| if on { c.abs().max(f64::MIN_POSITIVE) } else { 0.0 })
            .collect();
        // Radiating elements keep a tiny amplitude, so the limits never see an
        // all-zero excitation (a vanished one becomes uniform); the limits
        // were validated by `check`, so `apply` cannot fail.
        let amps = options.limits.apply(&amps).unwrap_or_default();
        w.iter()
            .zip(&amps)
            .enumerate()
            .map(|(n, (c, a))| {
                if !self.geo.on[n] {
                    return Complex::ZERO;
                }
                let s = self.geo.steer[n];
                let phase = quantise_phase(c.arg() + s, self.geo.phase_bits) - s;
                Complex::from_polar(*a, phase)
            })
            .collect()
    }
}

/// How the coverage level is set during the iterations.
#[derive(Debug, Clone, Copy)]
enum Mode {
    /// Fixed level (`None` without coverage-relative bounds).
    Fixed(Option<f64>),
    /// Refitted at each iteration to the coverage gains.
    Adaptive,
}

/// The better of two iterates.
fn better(a: Iterate, b: Iterate, tolerance_db: f64) -> Iterate {
    if a.better_than(&b, tolerance_db) { a } else { b }
}

/// Settings of one synthesis run.
struct Run<'a> {
    bounds: &'a [PatternBound],
    options: &'a SynthesisOptions,
}

/// An evaluated iterate with what the next step needs.
struct Evaluated {
    it: Iterate,
    field: Vec<Complex>,
    power: f64,
    target: Option<f64>,
}

/// One iterate of the solver.
#[derive(Clone)]
struct Iterate {
    w: Vec<Complex>,
    level: Option<f64>,
    violation: f64,
    gains: Vec<f64>,
}

impl Iterate {
    /// Both meet the bounds: the higher level wins; otherwise the smaller violation.
    fn better_than(&self, other: &Self, tolerance_db: f64) -> bool {
        if self.violation <= tolerance_db && other.violation <= tolerance_db {
            self.level.unwrap_or(0.0) > other.level.unwrap_or(0.0)
                || (self.level == other.level && self.violation < other.violation)
        } else {
            self.violation < other.violation
        }
    }
}

/// Median of `g − o_min` over the coverage samples (relative lower bounds).
fn median_coverage(bounds: &[PatternBound], gains: &[f64]) -> f64 {
    let mut v: Vec<f64> = bounds
        .iter()
        .zip(gains)
        .filter_map(|(b, g)| match b.min {
            Some(Level::Coverage(o)) => Some(g - o),
            _ => None,
        })
        .collect();
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

/// Whether some bound is a coverage-relative lower bound.
fn has_coverage_min(bounds: &[PatternBound]) -> bool {
    bounds.iter().any(|b| matches!(b.min, Some(Level::Coverage(_))))
}

/// Coverage levels compatible with the coverage samples (those with a
/// relative lower bound): at most `g − o_min` (the weakest sample), at least
/// `g − o_max` for those that also have a relative ripple bound
/// (`−∞` without one). Other samples never drive the level.
fn level_bracket(bounds: &[PatternBound], gains: &[f64]) -> (f64, f64) {
    bounds
        .iter()
        .zip(gains)
        .fold((f64::NEG_INFINITY, f64::INFINITY), |(lo, hi), (b, g)| {
            match (b.min, b.max) {
                (Some(Level::Coverage(min)), Some(Level::Coverage(max))) => (lo.max(g - max), hi.min(g - min)),
                (Some(Level::Coverage(min)), _) => (lo, hi.min(g - min)),
                _ => (lo, hi),
            }
        })
}

/// Largest violation (dB) of the bounds at coverage level `level`.
fn violation_at(bounds: &[PatternBound], gains: &[f64], level: Option<f64>) -> f64 {
    let c = level.unwrap_or(0.0);
    bounds.iter().zip(gains).fold(0.0f64, |acc, (b, g)| {
        let lo = db(b.min, c, f64::NEG_INFINITY);
        let hi = db(b.max, c, f64::INFINITY);
        acc.max(lo - g).max(g - hi)
    })
}

/// A level in dB at coverage level `c`, or `default` when absent.
fn db(l: Option<Level>, c: f64, default: f64) -> f64 {
    match l {
        None => default,
        Some(Level::Dbi(v)) => v,
        Some(Level::Coverage(o)) => c + o,
    }
}

/// A level as a linear gain at coverage level `c`, or `default` when absent.
fn linear(l: Option<Level>, c: f64, default: f64) -> f64 {
    match l {
        None => default,
        some => 10f64.powf(db(some, c, 0.0) / 10.0),
    }
}

fn null_steer(geo: &Geometry, dirs: &[[f64; 3]], w0: &[Complex]) -> Result<(Vec<f64>, Vec<f64>)> {
    if dirs.is_empty() {
        return Err(invalid("nulls", "give at least one direction"));
    }
    if dirs.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid("nulls", "directions must be finite"));
    }
    let radiating = geo.on.iter().filter(|o| **o).count();
    if dirs.len() >= radiating {
        return Err(invalid(
            "nulls",
            format!("{} nulls need more than {} radiating elements", dirs.len(), radiating),
        ));
    }
    let c: Vec<Vec<Complex>> = dirs
        .iter()
        .map(|d| {
            (0..w0.len())
                .map(|n| if geo.on[n] { geo.term(n, d) } else { Complex::ZERO })
                .collect()
        })
        .collect();
    // Gram matrix C Cᴴ and right-hand side C w₀.
    let k = c.len();
    let gram: Vec<Vec<Complex>> = (0..k).map(|i| (0..k).map(|j| dot(&c[i], &c[j])).collect()).collect();
    let rhs: Vec<Complex> = c
        .iter()
        .map(|ci| ci.iter().zip(w0).fold(Complex::ZERO, |s, (a, b)| s + *a * *b))
        .collect();
    let y = solve(gram, rhs)
        .ok_or_else(|| invalid("nulls", "the nulls cannot be met together (duplicate or too close)"))?;
    let w: Vec<Complex> = (0..w0.len())
        .map(|n| {
            if !geo.on[n] {
                return Complex::ZERO;
            }
            let corr = c
                .iter()
                .zip(&y)
                .fold(Complex::ZERO, |s, (ci, yi)| s + ci[n].conj() * *yi);
            w0[n] - corr
        })
        .collect();
    let max = w.iter().map(|c| c.abs()).fold(0.0, f64::max);
    if max <= 1e-12 * w0.iter().map(|c| c.abs()).fold(0.0, f64::max) {
        return Err(AntennaError::InvalidParameter {
            name: "nulls",
            reason: "the nulls remove the whole excitation".into(),
        });
    }
    Ok((
        w.iter().map(|c| c.abs() / max).collect(),
        w.iter().map(|c| c.arg()).collect(),
    ))
}

/// `Σ a_n conj(b_n)`.
fn dot(a: &[Complex], b: &[Complex]) -> Complex {
    a.iter().zip(b).fold(Complex::ZERO, |s, (x, y)| s + *x * y.conj())
}

/// Solves `A x = b` (Gaussian elimination, partial pivoting); `None` if singular.
fn solve(mut a: Vec<Vec<Complex>>, mut b: Vec<Complex>) -> Option<Vec<Complex>> {
    let n = b.len();
    let scale = a.iter().flatten().map(|c| c.abs()).fold(0.0, f64::max);
    for col in 0..n {
        let piv = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[piv][col].abs() <= 1e-10 * scale {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let inv = reciprocal(a[col][col]);
        for row in col + 1..n {
            let f = a[row][col] * inv;
            let pivot_row = a[col].clone();
            for (x, v) in a[row].iter_mut().zip(&pivot_row).skip(col) {
                *x -= f * *v;
            }
            let v = b[col];
            b[row] -= f * v;
        }
    }
    let mut x = vec![Complex::ZERO; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in row + 1..n {
            s -= a[row][k] * x[k];
        }
        x[row] = s * reciprocal(a[row][row]);
    }
    Some(x)
}

/// `1 / z` (non-zero `z`).
fn reciprocal(z: Complex) -> Complex {
    z.conj().scale(1.0 / z.norm_sqr())
}

fn norm_sqr(v: &[Complex]) -> f64 {
    v.iter().map(|c| c.norm_sqr()).sum::<f64>()
}

/// `Σ |w_n|²` over radiating elements (never 0).
fn power(w: &[Complex], on: &[bool]) -> f64 {
    w.iter()
        .zip(on)
        .filter(|(_, o)| **o)
        .map(|(c, _)| c.norm_sqr())
        .sum::<f64>()
        .max(MIN_GAIN_LIN)
}

/// `(0..n).map(f)`, in parallel with the `parallel` feature.
#[cfg(feature = "parallel")]
fn rows<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    use rayon::prelude::*;
    (0..n).into_par_iter().map(f).collect()
}

/// `(0..n).map(f)`, in parallel with the `parallel` feature.
#[cfg(not(feature = "parallel"))]
fn rows<T>(n: usize, f: impl Fn(usize) -> T) -> Vec<T> {
    (0..n).map(f).collect()
}

#[cfg(test)]
mod tests;
