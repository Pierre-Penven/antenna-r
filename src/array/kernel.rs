//! Shift-invariant gain kernel of a planar phased array: fast evaluation of
//! steered (natural or shaped) beams.
//!
//! For a planar array (`z = 0`) without phase quantisation, the per-element
//! weights `w_n = a_n e^{jφ_n}` (taper, excitation offsets, errors, failures)
//! do not depend on the steering, which only adds the linear phase
//! `−s (x_n u₀ + y_n v₀)` (`s = k₀` for phase shifters, `s = k` for true time
//! delay). The array factor at wavenumber `k` is then
//!
//! `AF(u, v) = Σ w_n e^{j k (x_n (u − σu₀) + y_n (v − σv₀))}`, `σ = s / k`,
//!
//! a function of the offset from the (squinted) beam only. The kernel tabulates
//! `|AF|² / Σa²` over the offsets once — as two 1D tables when the elements form
//! a rectangular grid with weights that factor by row and column, else as one
//! 2D table — so a gain costs a few table lookups instead of a sum over the
//! elements.
//!
//! ```
//! use antenna_r::{Direction, Isotropic};
//! use antenna_r::array::{PhasedArray, SteeringMode, lattice::rectangular};
//! let f = 20e9;
//! let d = rf_r::fspl::wavelength_m(f) / 2.0;
//! let mut a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &rectangular(20, 20, d, d).unwrap(), f).unwrap();
//! let kernel = a.kernel(f, 1e-4).unwrap();
//! assert!(kernel.is_separable());
//! // A beam steered to (u0, v0) = (0.2, -0.1): the kernel gives the array gain anywhere.
//! a.steer(Direction::from_uv(0.2, -0.1).unwrap(), SteeringMode::Phase);
//! let dir = Direction::from_uv(0.23, -0.08).unwrap();
//! let exact = a.array_gain_db(dir, f);
//! assert!((kernel.array_gain_db(0.23, -0.08, 0.2, -0.1) - exact).abs() < 1e-3);
//! ```

use crate::array::{PhasedArray, SteeringMode, wavenumber};
use crate::error::{Result, invalid, require_positive};
use crate::math::Complex;

/// Largest 2D table accepted (samples).
pub const MAX_TABLE_SAMPLES: usize = 64_000_000;

/// Relative tolerance on element coordinates and weights when detecting a
/// separable (row × column) array.
const SEPARABLE_TOL: f64 = 1e-9;

/// Samples of `f(lo + i·step)` with linear interpolation.
#[derive(Debug, Clone, PartialEq)]
struct Table1d {
    lo: f64,
    step: f64,
    values: Vec<f64>,
}

impl Table1d {
    #[inline]
    fn at(&self, x: f64) -> f64 {
        let p = ((x - self.lo) / self.step).clamp(0.0, (self.values.len() - 1) as f64);
        let i = (p as usize).min(self.values.len() - 2);
        let f = p - i as f64;
        self.values[i] + f * (self.values[i + 1] - self.values[i])
    }
}

/// Samples of `f(lo + i·step, lo + j·step)` (row `j`), bilinear interpolation.
#[derive(Debug, Clone, PartialEq)]
struct Table2d {
    lo: f64,
    step: f64,
    count: usize,
    values: Vec<f32>,
}

impl Table2d {
    #[inline]
    fn at(&self, x: f64, y: f64) -> f64 {
        let last = (self.count - 1) as f64;
        let px = ((x - self.lo) / self.step).clamp(0.0, last);
        let py = ((y - self.lo) / self.step).clamp(0.0, last);
        let (i, j) = ((px as usize).min(self.count - 2), (py as usize).min(self.count - 2));
        let (fx, fy) = (px - i as f64, py - j as f64);
        let v = |i: usize, j: usize| f64::from(self.values[j * self.count + i]);
        let top = v(i, j) + fx * (v(i + 1, j) - v(i, j));
        let bottom = v(i, j + 1) + fx * (v(i + 1, j + 1) - v(i, j + 1));
        top + fy * (bottom - top)
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Shape {
    /// `|AF|² = |AFx(Δu)|² · |AFy(Δv)|²`.
    Separable { x: Table1d, y: Table1d },
    /// `|AF(Δu, Δv)|²`.
    Grid(Table2d),
}

/// The shift-invariant array gain of a planar phased array at one frequency
/// (see the module documentation). The element pattern is not included.
#[derive(Debug, Clone, PartialEq)]
pub struct ArrayKernel {
    shape: Shape,
    /// `Σ a²`.
    power: f64,
    /// Beam shift per unit of steering: `k₀ / k` (phase) or 1 (true time delay).
    squint: f64,
}

impl ArrayKernel {
    /// Whether the kernel factors into two 1D tables (rectangular grid,
    /// row × column weights).
    #[must_use]
    pub fn is_separable(&self) -> bool {
        matches!(self.shape, Shape::Separable { .. })
    }

    /// Array gain `|AF|² / Σa²` (linear) towards direction cosines `(u, v)`
    /// of a beam steered to `(u0, v0)`.
    #[inline]
    #[must_use]
    pub fn array_gain(&self, u: f64, v: f64, u0: f64, v0: f64) -> f64 {
        let (du, dv) = (u - self.squint * u0, v - self.squint * v0);
        let af2 = match &self.shape {
            Shape::Separable { x, y } => x.at(du) * y.at(dv),
            Shape::Grid(t) => t.at(du, dv),
        };
        af2 / self.power
    }

    /// [`array_gain`](Self::array_gain) in dB (`−∞` in an exact null).
    #[inline]
    #[must_use]
    pub fn array_gain_db(&self, u: f64, v: f64, u0: f64, v0: f64) -> f64 {
        rf_r::db::to_db(self.array_gain(u, v, u0, v0))
    }
}

impl<E> PhasedArray<E> {
    /// The shift-invariant gain kernel of this array at `freq_hz`, tabulated
    /// every `step_uv` in direction-cosine offset (1D tables for a separable
    /// array, else a 2D table; choose a step well below the beamwidth, e.g.
    /// `λ / (40 L)` for an aperture of size `L`).
    ///
    /// # Errors
    /// [`crate::AntennaError::InvalidParameter`] for a non-positive frequency
    /// or step, a non-planar array (`z ≠ 0`), quantised phases (the steering
    /// phase then changes the beam shape), or a 2D table above
    /// [`MAX_TABLE_SAMPLES`].
    pub fn kernel(&self, freq_hz: f64, step_uv: f64) -> Result<ArrayKernel> {
        require_positive("freq_hz", freq_hz)?;
        require_positive("step_uv", step_uv)?;
        let g = &*self.geometry;
        if g.z.iter().any(|&z| z != 0.0) {
            return Err(invalid("positions", "the kernel needs a planar array (all z = 0)"));
        }
        if self.phase_bits.is_some() {
            return Err(invalid(
                "phase_bits",
                "quantised phases make the beam shape depend on the steering",
            ));
        }
        let k = wavenumber(freq_hz);
        let squint = match self.mode {
            SteeringMode::Phase => wavenumber(self.design_freq_hz) / k,
            SteeringMode::TrueTimeDelay => 1.0,
        };
        // Static weights: amplitude and the phase that does not come from the steering.
        let weights: Vec<Complex> = (0..self.amplitude.len())
            .map(|n| {
                let phase = self.excitation_phase_rad[n] + self.phase_errors_rad[n];
                Complex::new(self.amplitude[n] * phase.cos(), self.amplitude[n] * phase.sin())
            })
            .collect();
        // Offsets reach 1 + squint (|u| ≤ 1 and |σu₀| ≤ σ).
        let reach = 1.0 + squint;
        let shape = match separable(&g.x, &g.y, &weights) {
            Some((xs, a, ys, b)) => Shape::Separable {
                x: table_1d(&xs, &a, k, reach, step_uv),
                y: table_1d(&ys, &b, k, reach, step_uv),
            },
            None => Shape::Grid(table_2d(&g.x, &g.y, &weights, k, reach, step_uv)?),
        };
        Ok(ArrayKernel {
            shape,
            power: self.power,
            squint,
        })
    }
}

/// Distinct sorted values (within the tolerance) and the index of each input.
fn distinct(values: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let scale = values.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(f64::MIN_POSITIVE);
    let mut sorted: Vec<f64> = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut unique: Vec<f64> = Vec::new();
    for v in sorted {
        if unique.last().is_none_or(|&u| (v - u).abs() > SEPARABLE_TOL * scale) {
            unique.push(v);
        }
    }
    let index = values
        .iter()
        .map(|v| {
            unique
                .iter()
                .position(|u| (v - u).abs() <= SEPARABLE_TOL * scale)
                .unwrap_or(0)
        })
        .collect();
    (unique, index)
}

/// Column positions and weights, row positions and weights when the elements
/// fill a rectangular grid exactly once and `w(i, j) = a_i b_j`.
#[allow(clippy::type_complexity)]
fn separable(x: &[f64], y: &[f64], w: &[Complex]) -> Option<(Vec<f64>, Vec<Complex>, Vec<f64>, Vec<Complex>)> {
    let (xs, xi) = distinct(x);
    let (ys, yi) = distinct(y);
    if xs.len() * ys.len() != w.len() {
        return None;
    }
    let mut grid: Vec<Option<Complex>> = vec![None; w.len()];
    for n in 0..w.len() {
        let cell = &mut grid[yi[n] * xs.len() + xi[n]];
        if cell.is_some() {
            return None;
        }
        *cell = Some(w[n]);
    }
    let cell = |i: usize, j: usize| grid[j * xs.len() + i].unwrap_or_default();
    // Reference: the strongest element; rank one ⇔ w(i,j)·w(r) = w(i,jr)·w(ir,j).
    let r = (0..w.len()).max_by(|&a, &b| w[a].norm_sqr().total_cmp(&w[b].norm_sqr()))?;
    let (ir, jr) = (xi[r], yi[r]);
    let wr = w[r];
    let peak = wr.norm_sqr();
    for j in 0..ys.len() {
        for i in 0..xs.len() {
            let lhs = cell(i, j) * wr;
            let rhs = cell(i, jr) * cell(ir, j);
            if (lhs - rhs).norm_sqr() > (SEPARABLE_TOL * peak).powi(2) {
                return None;
            }
        }
    }
    let a: Vec<Complex> = (0..xs.len()).map(|i| cell(i, jr)).collect();
    let inv = Complex::new(wr.re / peak, -wr.im / peak);
    let b: Vec<Complex> = (0..ys.len()).map(|j| cell(ir, j) * inv).collect();
    Some((xs, a, ys, b))
}

/// `|Σ w_n e^{j k p_n Δ}|²` for `Δ` from `−reach` to `reach` every `step`,
/// by phasor recurrence (one complex multiply per element and sample).
fn row(positions: &[f64], weights: &[Complex], k: f64, lo: f64, step: f64, count: usize, out: &mut [f64]) {
    let mut phasors: Vec<Complex> = positions
        .iter()
        .zip(weights)
        .map(|(p, w)| {
            let a = k * p * lo;
            *w * Complex::new(a.cos(), a.sin())
        })
        .collect();
    let turns: Vec<Complex> = positions
        .iter()
        .map(|p| {
            let a = k * p * step;
            Complex::new(a.cos(), a.sin())
        })
        .collect();
    for (i, o) in out.iter_mut().enumerate().take(count) {
        let sum = phasors.iter().fold(Complex::default(), |s, &p| s + p);
        *o = sum.norm_sqr();
        // Re-anchor now and then so rounding does not accumulate.
        if i % 1024 == 1023 {
            let base = lo + (i + 1) as f64 * step;
            for ((ph, p), w) in phasors.iter_mut().zip(positions).zip(weights) {
                let a = k * p * base;
                *ph = *w * Complex::new(a.cos(), a.sin());
            }
        } else {
            for (ph, t) in phasors.iter_mut().zip(&turns) {
                *ph = *ph * *t;
            }
        }
    }
}

fn samples(reach: f64, step: f64) -> usize {
    (2.0 * reach / step).ceil() as usize + 2
}

fn table_1d(positions: &[f64], weights: &[Complex], k: f64, reach: f64, step: f64) -> Table1d {
    let count = samples(reach, step);
    let mut values = vec![0.0; count];
    row(positions, weights, k, -reach, step, count, &mut values);
    Table1d {
        lo: -reach,
        step,
        values,
    }
}

fn table_2d(x: &[f64], y: &[f64], weights: &[Complex], k: f64, reach: f64, step: f64) -> Result<Table2d> {
    let count = samples(reach, step);
    if count.saturating_mul(count) > MAX_TABLE_SAMPLES {
        return Err(invalid(
            "step_uv",
            format!("a 2D kernel table of {count}² samples is too large; use a larger step"),
        ));
    }
    let lo = -reach;
    // One row per Δv: weights carry the e^{j k y Δv} factor, then the Δu recurrence.
    let fill = |j: usize, out: &mut [f32]| {
        let dv = lo + j as f64 * step;
        let w: Vec<Complex> = weights
            .iter()
            .zip(y)
            .map(|(w, yy)| {
                let a = k * yy * dv;
                *w * Complex::new(a.cos(), a.sin())
            })
            .collect();
        let mut line = vec![0.0; count];
        row(x, &w, k, lo, step, count, &mut line);
        for (o, v) in out.iter_mut().zip(line) {
            *o = v as f32;
        }
    };
    let mut values = vec![0.0f32; count * count];
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        values
            .par_chunks_mut(count)
            .enumerate()
            .for_each(|(j, out)| fill(j, out));
    }
    #[cfg(not(feature = "parallel"))]
    for (j, out) in values.chunks_mut(count).enumerate() {
        fill(j, out);
    }
    Ok(Table2d {
        lo,
        step,
        count,
        values,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::array::lattice::{rectangular, triangular};
    use crate::{Direction, Isotropic};

    const F: f64 = 20e9;

    fn lambda() -> f64 {
        rf_r::fspl::wavelength_m(F)
    }

    fn square(n: usize) -> PhasedArray<Isotropic> {
        let d = lambda() / 2.0;
        PhasedArray::new(Isotropic::new(0.0).unwrap(), &rectangular(n, n, d, d).unwrap(), F).unwrap()
    }

    /// Largest |kernel − exact| (dB) over directions within `range_db` of the beam peak.
    fn max_error(a: &PhasedArray<Isotropic>, k: &ArrayKernel, freq: f64, u0: f64, v0: f64, range_db: f64) -> f64 {
        let peak = (0..=40)
            .flat_map(|i| (0..=40).map(move |j| (i, j)))
            .map(|(i, j)| (u0 - 0.2 + 0.01 * f64::from(i), v0 - 0.2 + 0.01 * f64::from(j)))
            .filter_map(|(u, v)| Direction::from_uv(u, v).ok())
            .map(|d| a.array_gain_db(d, freq))
            .fold(f64::NEG_INFINITY, f64::max);
        let mut worst = 0.0f64;
        for i in 0..60 {
            for j in 0..60 {
                let (u, v) = (-0.9 + 0.03 * f64::from(i), -0.9 + 0.03 * f64::from(j));
                let Ok(d) = Direction::from_uv(u, v) else { continue };
                let exact = a.array_gain_db(d, freq);
                if exact > peak - range_db {
                    worst = worst.max((k.array_gain_db(u, v, u0, v0) - exact).abs());
                }
            }
        }
        worst
    }

    #[test]
    fn a_uniform_square_array_is_separable_and_exact() {
        let mut a = square(20);
        let k = a.kernel(F, 5e-5).unwrap();
        assert!(k.is_separable());
        for (u0, v0) in [(0.0, 0.0), (0.2, -0.1), (-0.5, 0.4)] {
            a.steer(Direction::from_uv(u0, v0).unwrap(), SteeringMode::Phase);
            assert!(max_error(&a, &k, F, u0, v0, 50.0) < 0.05, "({u0}, {v0})");
        }
        // Linear gain at the beam: N for isotropic elements.
        assert!((k.array_gain(0.2, -0.1, 0.2, -0.1) - 400.0).abs() < 1e-6);
    }

    #[test]
    fn row_and_column_tapers_stay_separable() {
        let n = 12;
        let a_row: Vec<f64> = (0..n).map(|i| 1.0 - 0.5 * ((i as f64 - 5.5) / 6.0).powi(2)).collect();
        let amps: Vec<f64> = (0..n * n).map(|m| a_row[m % n] * a_row[m / n]).collect();
        let mut a = square(n).with_amplitudes(&amps).unwrap();
        let k = a.kernel(F, 5e-5).unwrap();
        assert!(k.is_separable());
        a.steer(Direction::from_uv(0.3, 0.1).unwrap(), SteeringMode::Phase);
        assert!(max_error(&a, &k, F, 0.3, 0.1, 50.0) < 0.05);
    }

    #[test]
    fn phase_steering_squints_off_frequency_and_true_time_delay_does_not() {
        let mut a = square(16);
        a.steer(Direction::from_uv(0.4, 0.0).unwrap(), SteeringMode::Phase);
        let f = 1.1 * F;
        let k = a.kernel(f, 5e-5).unwrap();
        assert!(max_error(&a, &k, f, 0.4, 0.0, 40.0) < 0.05);
        // The beam peaks at u = (f0 / f) u0.
        let squinted = 0.4 / 1.1;
        assert!(k.array_gain(squinted, 0.0, 0.4, 0.0) > k.array_gain(0.4, 0.0, 0.4, 0.0));
        a.steer(Direction::from_uv(0.4, 0.0).unwrap(), SteeringMode::TrueTimeDelay);
        let t = a.kernel(f, 5e-5).unwrap();
        assert!(max_error(&a, &t, f, 0.4, 0.0, 40.0) < 0.05);
        assert!((t.array_gain(0.4, 0.0, 0.4, 0.0) - 256.0).abs() < 1e-6);
    }

    #[test]
    fn shaped_and_irregular_arrays_use_a_2d_table() {
        // Triangular lattice: not a rectangular grid.
        let d = lambda() / 2.0;
        let pos = triangular(10, 10, d, d * 0.866).unwrap();
        let mut a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &pos, F).unwrap();
        let k = a.kernel(F, 1e-3).unwrap();
        assert!(!k.is_separable());
        a.steer(Direction::from_uv(-0.2, 0.3).unwrap(), SteeringMode::Phase);
        assert!(max_error(&a, &k, F, -0.2, 0.3, 30.0) < 0.1);
        // Square grid with non-separable excitation (a shaped beam): 2D, still shift invariant.
        let n = 8;
        let phases: Vec<f64> = (0..n * n).map(|m| 0.3 * ((m % n) * (m / n)) as f64).collect();
        let mut s = square(n).with_excitation(&vec![1.0; n * n], &phases).unwrap();
        let ks = s.kernel(F, 1e-3).unwrap();
        assert!(!ks.is_separable());
        s.steer(Direction::from_uv(0.1, 0.1).unwrap(), SteeringMode::Phase);
        assert!(max_error(&s, &ks, F, 0.1, 0.1, 30.0) < 0.1);
        // Failed element: breaks the row × column product of a uniform grid.
        let f = square(n).with_failed_elements(&[9]).unwrap();
        assert!(!f.kernel(F, 2e-3).unwrap().is_separable());
    }

    #[test]
    fn rejects_what_is_not_shift_invariant_or_too_large() {
        let a = square(4);
        assert!(a.kernel(0.0, 1e-3).is_err());
        assert!(a.kernel(F, 0.0).is_err());
        assert!(a.clone().with_phase_bits(Some(4)).unwrap().kernel(F, 1e-3).is_err());
        let lifted = vec![
            geometry_r::math::Vec3::new(0.0, 0.0, 0.0),
            geometry_r::math::Vec3::new(0.01, 0.0, 0.005),
        ];
        let z = PhasedArray::new(Isotropic::new(0.0).unwrap(), &lifted, F).unwrap();
        assert!(z.kernel(F, 1e-3).is_err());
        let d = lambda() / 2.0;
        let tri = PhasedArray::new(Isotropic::new(0.0).unwrap(), &triangular(4, 4, d, d).unwrap(), F).unwrap();
        assert!(tri.kernel(F, 1e-5).is_err());
    }

    #[test]
    fn tables_clamp_outside_their_range_and_detect_duplicates() {
        let t = Table1d {
            lo: 0.0,
            step: 1.0,
            values: vec![1.0, 3.0],
        };
        assert!((t.at(-5.0) - 1.0).abs() < 1e-12);
        assert!((t.at(0.5) - 2.0).abs() < 1e-12);
        assert!((t.at(9.0) - 3.0).abs() < 1e-12);
        let g = Table2d {
            lo: 0.0,
            step: 1.0,
            count: 2,
            values: vec![0.0, 1.0, 2.0, 3.0],
        };
        assert!((g.at(0.5, 0.5) - 1.5).abs() < 1e-6);
        assert!((g.at(-1.0, 9.0) - 2.0).abs() < 1e-6);
        // Two elements at the same place: not a grid.
        let w = vec![Complex::new(1.0, 0.0); 2];
        assert!(separable(&[0.0, 0.0], &[0.0, 0.0], &w).is_none());
        // A full grid with a duplicate and a hole.
        let w4 = vec![Complex::new(1.0, 0.0); 4];
        assert!(separable(&[0.0, 1.0, 0.0, 0.0], &[0.0, 0.0, 1.0, 1.0], &w4).is_none());
        // Long rows re-anchor the recurrence without drift.
        let mut out = vec![0.0; 3000];
        row(
            &[0.0, 0.01],
            &[Complex::new(1.0, 0.0); 2],
            100.0,
            -2.0,
            1e-3,
            3000,
            &mut out,
        );
        let exact = |dl: f64| {
            let a = 100.0 * 0.01 * dl;
            (1.0 + a.cos()).powi(2) + a.sin().powi(2)
        };
        assert!((out[2999] - exact(-2.0 + 2999.0 * 1e-3)).abs() < 1e-9);
    }
}
