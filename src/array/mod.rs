//! Phased arrays: element pattern × array factor.
//!
//! ```text
//! G(û, f) = G_e(û) · |AF(û, f)|² / Σ a_n²
//! AF(û, f) = Σ a_n exp(j (k r_n · û + ψ_n))
//! ```
//!
//! `G_e` is the embedded (in-array) element gain, `a_n` the amplitudes,
//! `r_n` the element positions and `ψ_n` the excitation phases. The
//! normalisation by `Σ a_n²` gives the array gain relative to the power fed
//! to the radiating elements and is exact for well-sampled lattices without
//! mutual coupling (ADR 0004); failed elements radiate nothing and are left
//! out of the sum. Use [`crate::analysis::directivity_dbi`] for an exact
//! numerical directivity.
//!
//! Steering: [`SteeringMode::Phase`] sets `ψ_n = −k0 r_n · û0` at the
//! design frequency, so the beam squints at other frequencies;
//! [`SteeringMode::TrueTimeDelay`] uses `−k r_n · û0` at every frequency.

pub mod excitation;
pub mod grating;
pub mod lattice;
pub mod synthesis;
pub mod taper;

use std::f64::consts::TAU;
use std::sync::Arc;

use geometry_r::math::Vec3;

use crate::direction::Direction;
use crate::error::{AntennaError, Result, invalid, require_level_db, require_positive};
use crate::math::Complex;
use crate::pattern::{DEFAULT_FLOOR_DBI, GainPattern};

/// How the steering phases depend on frequency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SteeringMode {
    /// Phase shifters set at the design frequency (beam squint).
    #[default]
    Phase,
    /// True time delay (no squint).
    TrueTimeDelay,
}

/// Element positions, shared between the beams of a multi-beam array.
#[derive(Debug, PartialEq)]
struct Geometry {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

/// A phased array of identical elements.
#[derive(Debug, Clone, PartialEq)]
pub struct PhasedArray<E> {
    element: E,
    geometry: Arc<Geometry>,
    design_freq_hz: f64,
    taper: Vec<f64>,
    excitation_phase_rad: Vec<f64>,
    failed: Vec<bool>,
    amp_errors_db: Vec<f64>,
    phase_errors_rad: Vec<f64>,
    phase_bits: Option<u8>,
    steering: Direction,
    mode: SteeringMode,
    floor_dbi: f64,
    // Derived state, rebuilt by `rebuild`.
    amplitude: Vec<f64>,
    projection: Vec<f64>,
    phase: Vec<f64>,
    power: f64,
}

/// Wavenumber `2π f / c` (rad/m).
fn wavenumber(freq_hz: f64) -> f64 {
    TAU * freq_hz / rf_r::constants::SPEED_OF_LIGHT
}

impl<E: GainPattern> PhasedArray<E> {
    /// Array of `element` at `positions` (metres, antenna frame), uniformly
    /// excited and steered to broadside, for `design_freq_hz`.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, Isotropic};
    /// use antenna_r::array::{PhasedArray, lattice::linear};
    /// let lambda = rf_r::fspl::wavelength_m(10e9);
    /// let pos = linear(8, lambda / 2.0).unwrap();
    /// let a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &pos, 10e9).unwrap();
    /// assert!((a.gain_dbi(Direction::BORESIGHT) - 10.0 * 8f64.log10()).abs() < 1e-9);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for no positions, a non-finite
    /// position or a non-positive frequency.
    pub fn new(element: E, positions: &[Vec3], design_freq_hz: f64) -> Result<Self> {
        if positions.is_empty() {
            return Err(invalid("positions", "an array needs at least one element"));
        }
        if positions.iter().any(|p| !p.is_finite()) {
            return Err(invalid("positions", "positions must be finite"));
        }
        let n = positions.len();
        let mut array = Self {
            element,
            geometry: Arc::new(Geometry {
                x: positions.iter().map(|p| p.x).collect(),
                y: positions.iter().map(|p| p.y).collect(),
                z: positions.iter().map(|p| p.z).collect(),
            }),
            design_freq_hz: require_positive("design_freq_hz", design_freq_hz)?,
            taper: vec![1.0; n],
            excitation_phase_rad: vec![0.0; n],
            failed: vec![false; n],
            amp_errors_db: vec![0.0; n],
            phase_errors_rad: vec![0.0; n],
            phase_bits: None,
            steering: Direction::BORESIGHT,
            mode: SteeringMode::Phase,
            floor_dbi: DEFAULT_FLOOR_DBI,
            amplitude: Vec::new(),
            projection: Vec::new(),
            phase: Vec::new(),
            power: 0.0,
        };
        array.rebuild()?;
        Ok(array)
    }

    /// Number of elements (failed ones included).
    #[must_use]
    pub fn len(&self) -> usize {
        self.taper.len()
    }

    /// Always `false`: an array has at least one element.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Design frequency (Hz).
    #[must_use]
    pub fn design_freq_hz(&self) -> f64 {
        self.design_freq_hz
    }

    /// Commanded beam direction and steering mode.
    #[must_use]
    pub fn steering(&self) -> (Direction, SteeringMode) {
        (self.steering, self.mode)
    }

    /// The element pattern.
    #[must_use]
    pub fn element(&self) -> &E {
        &self.element
    }

    /// Element amplitudes (taper), one per element, all ≥ 0, not all zero.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for a length mismatch, negative or
    /// non-finite values, or no radiating element.
    pub fn with_amplitudes(mut self, amplitudes: &[f64]) -> Result<Self> {
        self.check_len("amplitudes", amplitudes.len())?;
        if amplitudes.iter().any(|a| !(a.is_finite() && *a >= 0.0)) {
            return Err(invalid("amplitudes", "must be finite and ≥ 0"));
        }
        self.taper = amplitudes.to_vec();
        self.rebuild()?;
        Ok(self)
    }

    /// Per-element excitation: amplitudes (as [`with_amplitudes`](Self::with_amplitudes))
    /// and phase offsets added to the steering phase, so a beam shape
    /// designed at broadside is steered as a whole. With
    /// [`with_phase_bits`](Self::with_phase_bits) the *total* phase
    /// (steering + offset) is quantised. Element `n` radiates
    /// `a_n exp(j(φ_steer,n + φ_n))`.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, Isotropic};
    /// use antenna_r::array::{PhasedArray, lattice::linear};
    /// let lambda = rf_r::fspl::wavelength_m(10e9);
    /// let pos = linear(2, lambda / 2.0).unwrap();
    /// // Opposite phases: a null on broadside.
    /// let a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &pos, 10e9)
    ///     .unwrap()
    ///     .with_excitation(&[1.0, 1.0], &[0.0, std::f64::consts::PI])
    ///     .unwrap()
    ///     .with_floor_dbi(-200.0)
    ///     .unwrap();
    /// assert!(a.gain_dbi(Direction::BORESIGHT) < -100.0);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for length mismatches, negative or
    /// non-finite amplitudes, non-finite phases, or no radiating element.
    pub fn with_excitation(mut self, amplitudes: &[f64], phases_rad: &[f64]) -> Result<Self> {
        self.check_len("phases_rad", phases_rad.len())?;
        if phases_rad.iter().any(|p| !p.is_finite()) {
            return Err(invalid("phases_rad", "must be finite"));
        }
        self.excitation_phase_rad = phases_rad.to_vec();
        self.with_amplitudes(amplitudes)
    }

    /// Phase quantisation (bits), if any.
    #[must_use]
    pub fn phase_bits(&self) -> Option<u8> {
        self.phase_bits
    }

    /// Excitation amplitudes (before failures and errors) and phase offsets (rad).
    #[must_use]
    pub fn excitation(&self) -> (&[f64], &[f64]) {
        (&self.taper, &self.excitation_phase_rad)
    }

    /// Element positions (m, antenna frame).
    #[must_use]
    pub fn positions(&self) -> Vec<Vec3> {
        let g = &*self.geometry;
        (0..self.len()).map(|n| Vec3::new(g.x[n], g.y[n], g.z[n])).collect()
    }

    /// Failed-element flags, one per element.
    #[must_use]
    pub fn failed(&self) -> &[bool] {
        &self.failed
    }

    /// Ideal steering phase of each element at the design frequency (rad,
    /// `−k₀ r_n · ŝ`), before offsets and quantisation.
    #[must_use]
    pub fn steering_phases_rad(&self) -> Vec<f64> {
        let k0 = wavenumber(self.design_freq_hz);
        self.projection.iter().map(|p| -k0 * p).collect()
    }

    /// Radiated amplitudes actually used (failures and amplitude errors applied).
    #[must_use]
    pub fn radiated_amplitudes(&self) -> &[f64] {
        &self.amplitude
    }

    /// Power efficiency of one amplifier per element at a fixed peak drive:
    /// `Σa² / (N · max a²)` (1 for uniform excitation). The EIRP of a
    /// transmit array is lower than the uniform one by this factor on top of
    /// the directivity change; a receive array or a passive feed only sees
    /// the directivity. See [`taper::power_efficiency`].
    #[must_use]
    pub fn power_efficiency(&self) -> f64 {
        taper::power_efficiency(&self.amplitude)
    }

    /// Marks elements as failed (they radiate nothing).
    ///
    /// # Errors
    /// [`AntennaError::OutOfRange`] for an index ≥ `len()`;
    /// [`AntennaError::InvalidParameter`] if no element would radiate.
    pub fn with_failed_elements(mut self, indices: &[usize]) -> Result<Self> {
        let mut failed = vec![false; self.len()];
        for &i in indices {
            *failed
                .get_mut(i)
                .ok_or_else(|| AntennaError::OutOfRange(format!("element {i} of {}", self.len())))? = true;
        }
        self.failed = failed;
        self.rebuild()?;
        Ok(self)
    }

    /// Per-element amplitude (dB) and phase (rad) errors, e.g. drawn by the
    /// caller from a seeded random generator.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for length mismatches or
    /// non-finite values.
    pub fn with_errors(mut self, amplitude_db: &[f64], phase_rad: &[f64]) -> Result<Self> {
        self.check_len("amplitude_db", amplitude_db.len())?;
        self.check_len("phase_rad", phase_rad.len())?;
        if amplitude_db.iter().chain(phase_rad).any(|v| !v.is_finite()) {
            return Err(invalid("errors", "must be finite"));
        }
        self.amp_errors_db = amplitude_db.to_vec();
        self.phase_errors_rad = phase_rad.to_vec();
        self.rebuild()?;
        Ok(self)
    }

    /// Quantises the phase-shifter settings to `bits` bits (1–16), or removes
    /// the quantisation with `None`. In [`SteeringMode::Phase`] the total
    /// phase (steering + excitation offset) is quantised; with true time
    /// delay only the excitation offsets are phases.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for `bits` outside 1–16.
    pub fn with_phase_bits(mut self, bits: Option<u8>) -> Result<Self> {
        if let Some(b) = bits
            && !(1..=16).contains(&b)
        {
            return Err(invalid("bits", format!("must be in 1..=16, got {b}")));
        }
        self.phase_bits = bits;
        self.rebuild()?;
        Ok(self)
    }

    /// Lowest gain returned (nulls), in dBi.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] unless finite or `-inf`.
    pub fn with_floor_dbi(mut self, floor_dbi: f64) -> Result<Self> {
        self.floor_dbi = require_level_db("floor_dbi", floor_dbi)?;
        Ok(self)
    }

    /// Steers the beam towards `target` (antenna frame).
    pub fn steer(&mut self, target: Direction, mode: SteeringMode) {
        self.steering = target;
        self.mode = mode;
        self.update_phases();
    }

    /// Taper efficiency of the radiating elements.
    #[must_use]
    pub fn taper_efficiency(&self) -> f64 {
        taper::taper_efficiency(&self.amplitude)
    }

    /// Array factor `AF(û, f)` (complex, not normalised).
    #[must_use]
    pub fn array_factor(&self, dir: Direction, freq_hz: f64) -> Complex {
        let k = wavenumber(freq_hz);
        let ttd = if self.mode == SteeringMode::TrueTimeDelay {
            k
        } else {
            0.0
        };
        let (u, v, w) = (dir.u(), dir.v(), dir.w());
        let g = &*self.geometry;
        let mut re = 0.0;
        let mut im = 0.0;
        for n in 0..self.amplitude.len() {
            let path = g.x[n].mul_add(u, g.y[n].mul_add(v, g.z[n] * w));
            let phase = k.mul_add(path, ttd.mul_add(-self.projection[n], self.phase[n]));
            let (s, c) = phase.sin_cos();
            re += self.amplitude[n] * c;
            im += self.amplitude[n] * s;
        }
        Complex::new(re, im)
    }

    /// Array gain `10 log10(|AF|² / Σ a²)` in dB (element pattern excluded).
    #[must_use]
    pub fn array_gain_db(&self, dir: Direction, freq_hz: f64) -> f64 {
        rf_r::db::to_db(self.array_factor(dir, freq_hz).norm_sqr() / self.power)
    }

    /// Gain in dBi towards `dir` at `freq_hz`.
    #[must_use]
    pub fn gain_dbi_at(&self, dir: Direction, freq_hz: f64) -> f64 {
        if dir.w().is_nan() {
            return f64::NAN;
        }
        (self.element.gain_dbi(dir) + self.array_gain_db(dir, freq_hz)).max(self.floor_dbi)
    }

    fn check_len(&self, name: &'static str, len: usize) -> Result<()> {
        if len == self.len() {
            Ok(())
        } else {
            Err(invalid(name, format!("expected {} values, got {len}", self.len())))
        }
    }

    fn rebuild(&mut self) -> Result<()> {
        self.amplitude = self
            .taper
            .iter()
            .zip(&self.failed)
            .zip(&self.amp_errors_db)
            .map(|((a, f), e)| if *f { 0.0 } else { a * rf_r::db::amplitude_from_db(*e) })
            .collect();
        self.power = self.amplitude.iter().map(|a| a * a).sum();
        if self.power <= 0.0 {
            return Err(invalid("amplitudes", "no radiating element"));
        }
        self.update_phases();
        Ok(())
    }

    fn update_phases(&mut self) {
        let (u0, v0, w0) = (self.steering.u(), self.steering.v(), self.steering.w());
        let g = &*self.geometry;
        self.projection = (0..self.taper.len())
            .map(|n| g.x[n].mul_add(u0, g.y[n].mul_add(v0, g.z[n] * w0)))
            .collect();
        let k0 = wavenumber(self.design_freq_hz);
        let quantum = self.phase_bits.map(|b| TAU / f64::from(1u32 << b));
        self.phase = self
            .projection
            .iter()
            .zip(&self.excitation_phase_rad)
            .zip(&self.phase_errors_rad)
            .map(|((p, x), e)| {
                // Phase shifters set the total phase (steering + offset);
                // true time delay steers by delay, the offset stays a phase.
                let ideal = match self.mode {
                    SteeringMode::Phase => -k0 * p + x,
                    SteeringMode::TrueTimeDelay => *x,
                };
                quantum.map_or(ideal, |q| (ideal / q).round() * q) + e
            })
            .collect();
    }
}

impl<E: GainPattern> GainPattern for PhasedArray<E> {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_dbi_at(dir, self.design_freq_hz)
    }

    /// Upper bound of the gain: element peak + `10 log10((Σ a)² / Σ a²)`
    /// (all element contributions in phase).
    fn peak_gain_dbi(&self) -> f64 {
        let sum: f64 = self.amplitude.iter().sum();
        self.element.peak_gain_dbi() + rf_r::db::to_db(sum * sum / self.power)
    }
}

/// Several simultaneous beams sharing the element positions of one array.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiBeamArray<E> {
    beams: Vec<PhasedArray<E>>,
}

impl<E: GainPattern + Clone> MultiBeamArray<E> {
    /// One beam per steering of `base` (taper, failures and errors shared).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for an empty steering list.
    pub fn new(base: &PhasedArray<E>, steerings: &[(Direction, SteeringMode)]) -> Result<Self> {
        if steerings.is_empty() {
            return Err(invalid("steerings", "need at least one beam"));
        }
        let beams = steerings
            .iter()
            .map(|(d, m)| {
                let mut b = base.clone();
                b.steer(*d, *m);
                b
            })
            .collect();
        Ok(Self { beams })
    }

    /// Number of beams.
    #[must_use]
    pub fn beam_count(&self) -> usize {
        self.beams.len()
    }

    /// Beam `index`, or `None`.
    #[must_use]
    pub fn beam(&self, index: usize) -> Option<&PhasedArray<E>> {
        self.beams.get(index)
    }

    /// Re-steers beam `index`.
    ///
    /// # Errors
    /// [`AntennaError::OutOfRange`] for a missing beam.
    pub fn steer(&mut self, index: usize, target: Direction, mode: SteeringMode) -> Result<()> {
        let n = self.beams.len();
        self.beams
            .get_mut(index)
            .ok_or_else(|| AntennaError::OutOfRange(format!("beam {index} of {n}")))?
            .steer(target, mode);
        Ok(())
    }

    /// Gain of every beam towards `dir` at the design frequency.
    ///
    /// # Panics
    /// Panics if `out.len() != beam_count()`.
    pub fn gains_dbi(&self, dir: Direction, out: &mut [f64]) {
        assert_eq!(out.len(), self.beams.len(), "out must hold one value per beam");
        for (o, b) in out.iter_mut().zip(&self.beams) {
            *o = b.gain_dbi(dir);
        }
    }

    /// Whether all beams share the same element positions.
    #[must_use]
    pub fn shares_geometry(&self) -> bool {
        self.beams
            .windows(2)
            .all(|w| Arc::ptr_eq(&w[0].geometry, &w[1].geometry))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::mean_gain_db_piecewise;
    use crate::math::assert_close;
    use crate::pattern::{CosineQ, Isotropic};
    use std::f64::consts::PI;

    const F: f64 = 10e9;

    fn lambda() -> f64 {
        rf_r::fspl::wavelength_m(F)
    }

    fn iso() -> Isotropic {
        Isotropic::new(0.0).unwrap()
    }

    fn ula(n: usize) -> PhasedArray<Isotropic> {
        PhasedArray::new(iso(), &lattice::linear(n, lambda() / 2.0).unwrap(), F).unwrap()
    }

    #[test]
    fn broadside_gain_of_uniform_linear_array_is_ten_log_n() {
        for n in [1, 4, 16] {
            let a = ula(n);
            assert_close(a.gain_dbi(Direction::BORESIGHT), 10.0 * (n as f64).log10(), 1e-9);
            assert_close(a.peak_gain_dbi(), 10.0 * (n as f64).log10(), 1e-9);
            assert_eq!(a.len(), n);
            assert!(!a.is_empty());
        }
    }

    #[test]
    fn half_wave_linear_array_directivity_is_n() {
        // λ/2 isotropic ULA: exact directivity N at broadside.
        let a = ula(8);
        let mean = mean_gain_db_piecewise(&a, &[0.0, PI / 2.0, PI], 200, 360);
        assert_close(mean, 0.0, 1e-6);
    }

    #[test]
    fn phase_steering_moves_the_beam_and_squints_off_frequency() {
        let mut a = ula(32);
        let theta0 = 30f64.to_radians();
        a.steer(Direction::from_theta_phi(theta0, 0.0), SteeringMode::Phase);
        assert_close(
            a.gain_dbi(Direction::from_theta_phi(theta0, 0.0)),
            10.0 * 32f64.log10(),
            1e-9,
        );
        assert_eq!(a.steering().1, SteeringMode::Phase);
        // At 1.1 f0 the beam peaks at sin θ = (f0/f) sin θ0.
        let f = 1.1 * F;
        let squint = (theta0.sin() / 1.1).asin();
        let at_squint = a.gain_dbi_at(Direction::from_theta_phi(squint, 0.0), f);
        assert_close(at_squint, 10.0 * 32f64.log10(), 1e-9);
        assert!(a.gain_dbi_at(Direction::from_theta_phi(theta0, 0.0), f) < at_squint - 1.0);
    }

    #[test]
    fn true_time_delay_does_not_squint() {
        let mut a = ula(32);
        let target = Direction::from_theta_phi(30f64.to_radians(), 0.0);
        a.steer(target, SteeringMode::TrueTimeDelay);
        for f in [0.8 * F, F, 1.2 * F] {
            assert_close(a.gain_dbi_at(target, f), 10.0 * 32f64.log10(), 1e-9);
        }
    }

    #[test]
    fn element_pattern_multiplies_the_array_factor() {
        let e = CosineQ::new(1.0).unwrap();
        let a = PhasedArray::new(e, &lattice::linear(4, lambda() / 2.0).unwrap(), F).unwrap();
        let d = Direction::from_theta_phi(0.2, PI / 2.0); // across the array axis: AF = N
        assert_close(a.gain_dbi(d), e.gain_dbi(d) + 10.0 * 4f64.log10(), 1e-9);
        assert_eq!(a.element().peak_gain_dbi(), e.peak_gain_dbi());
    }

    #[test]
    fn taper_reduces_peak_by_its_efficiency() {
        let w = taper::chebyshev(16, 30.0).unwrap();
        let a = ula(16).with_amplitudes(&w).unwrap();
        let expected = 10.0 * (16.0 * taper::taper_efficiency(&w)).log10();
        assert_close(a.gain_dbi(Direction::BORESIGHT), expected, 1e-9);
        assert_close(a.taper_efficiency(), taper::taper_efficiency(&w), 1e-12);
    }

    #[test]
    fn failed_elements_and_errors() {
        let a = ula(4).with_failed_elements(&[0]).unwrap();
        assert_close(a.gain_dbi(Direction::BORESIGHT), 10.0 * 3f64.log10(), 1e-9);
        let e = ula(4).with_errors(&[0.0; 4], &[PI, 0.0, 0.0, 0.0]).unwrap();
        // One of four elements in anti-phase: |AF|² = 4, Σa² = 4 → 0 dB.
        assert_close(e.gain_dbi(Direction::BORESIGHT), 0.0, 1e-9);
        let amp = ula(2).with_errors(&[20.0 * 2f64.log10(), 0.0], &[0.0, 0.0]).unwrap();
        // Amplitudes 2 and 1: |AF|² = 9, Σa² = 5.
        assert_close(amp.gain_dbi(Direction::BORESIGHT), 10.0 * 1.8f64.log10(), 1e-9);
    }

    #[test]
    fn phase_quantisation_matches_bits() {
        let mut a = ula(8).with_phase_bits(Some(2)).unwrap();
        let target = Direction::from_theta_phi(10f64.to_radians(), 0.0);
        a.steer(target, SteeringMode::Phase);
        let q = PI / 2.0;
        assert!(a.phase.iter().all(|p| ((p / q) - (p / q).round()).abs() < 1e-9));
        assert!(a.gain_dbi(target) < 10.0 * 8f64.log10());
        let b = a.clone().with_phase_bits(None).unwrap();
        assert_close(b.gain_dbi(target), 10.0 * 8f64.log10(), 1e-9);
        assert!(ula(2).with_phase_bits(Some(0)).is_err());
        assert!(ula(2).with_phase_bits(Some(17)).is_err());
    }

    #[test]
    fn floor_and_nan_handling() {
        let a = ula(2).with_floor_dbi(-30.0).unwrap();
        // Two elements λ/2 apart: null at endfire (u = 1).
        assert_eq!(a.gain_dbi(Direction::from_theta_phi(PI / 2.0, 0.0)), -30.0);
        assert!(a.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
        assert!(ula(2).with_floor_dbi(f64::NAN).is_err());
    }

    #[test]
    fn construction_and_setters_validate() {
        assert!(matches!(
            PhasedArray::new(iso(), &[], F),
            Err(AntennaError::InvalidParameter { name: "positions", .. })
        ));
        assert!(matches!(
            PhasedArray::new(iso(), &[Vec3::new(f64::NAN, 0.0, 0.0)], F),
            Err(AntennaError::InvalidParameter { name: "positions", .. })
        ));
        assert!(PhasedArray::new(iso(), &[Vec3::ZERO], 0.0).is_err());
        assert!(matches!(
            ula(2).with_amplitudes(&[1.0]),
            Err(AntennaError::InvalidParameter { name: "amplitudes", .. })
        ));
        assert!(ula(2).with_amplitudes(&[1.0, -1.0]).is_err());
        assert!(ula(2).with_amplitudes(&[0.0, 0.0]).is_err());
        assert!(matches!(
            ula(2).with_failed_elements(&[2]),
            Err(AntennaError::OutOfRange(_))
        ));
        assert!(ula(2).with_failed_elements(&[0, 1]).is_err());
        assert!(ula(2).with_errors(&[0.0], &[0.0, 0.0]).is_err());
        assert!(ula(2).with_errors(&[0.0, 0.0], &[0.0]).is_err());
        assert!(ula(2).with_errors(&[0.0, f64::NAN], &[0.0, 0.0]).is_err());
        assert_eq!(ula(2).design_freq_hz(), F);
    }

    #[test]
    fn array_factor_of_one_element_is_its_amplitude() {
        let a = PhasedArray::new(iso(), &[Vec3::new(0.3, 0.1, 0.0)], F).unwrap();
        let af = a.array_factor(Direction::from_theta_phi(0.4, 1.0), F);
        assert_close(af.norm_sqr(), 1.0, 1e-12);
        assert_close(a.array_gain_db(Direction::BORESIGHT, F), 0.0, 1e-12);
    }

    #[test]
    fn multi_beam_shares_geometry_and_steers_each_beam() {
        let base = ula(16);
        let d1 = Direction::from_theta_phi(0.2, 0.0);
        let d2 = Direction::from_theta_phi(0.4, PI);
        let mut mb =
            MultiBeamArray::new(&base, &[(d1, SteeringMode::Phase), (d2, SteeringMode::TrueTimeDelay)]).unwrap();
        assert_eq!(mb.beam_count(), 2);
        assert!(mb.shares_geometry());
        let mut out = [0.0; 2];
        mb.gains_dbi(d1, &mut out);
        assert_close(out[0], 10.0 * 16f64.log10(), 1e-9);
        assert!(out[1] < out[0]);
        mb.steer(1, d1, SteeringMode::Phase).unwrap();
        assert_close(mb.beam(1).unwrap().gain_dbi(d1), 10.0 * 16f64.log10(), 1e-9);
        assert!(mb.beam(2).is_none());
        assert!(matches!(
            mb.steer(5, d1, SteeringMode::Phase),
            Err(AntennaError::OutOfRange(_))
        ));
        assert!(MultiBeamArray::new(&base, &[]).is_err());
    }

    #[test]
    #[should_panic(expected = "one value per beam")]
    fn multi_beam_gains_panic_on_bad_buffer() {
        let mb = MultiBeamArray::new(&ula(2), &[(Direction::BORESIGHT, SteeringMode::Phase)]).unwrap();
        mb.gains_dbi(Direction::BORESIGHT, &mut []);
    }

    #[test]
    fn excitation_offsets_add_to_the_steering_and_are_quantised_in_total() {
        let mut a = ula(4);
        a.steer(Direction::from_theta_phi(0.3, 0.0), SteeringMode::Phase);
        let offsets = [0.1, -0.2, 0.3, 0.0];
        let b = a.clone().with_excitation(&[1.0, 0.5, 1.0, 0.0], &offsets).unwrap();
        let (amps, phases) = b.excitation();
        assert_eq!(amps, &[1.0, 0.5, 1.0, 0.0]);
        assert_eq!(phases, &offsets);
        let steer = b.steering_phases_rad();
        for n in 0..4 {
            assert_close(b.phase[n], steer[n] + offsets[n], 1e-12);
        }
        assert_eq!(b.radiated_amplitudes(), &[1.0, 0.5, 1.0, 0.0]);
        assert_close(b.power_efficiency(), 2.25 / 4.0, 1e-12);
        // Quantisation applies to the total phase.
        let q = b.clone().with_phase_bits(Some(3)).unwrap();
        assert_eq!(q.phase_bits(), Some(3));
        let quantum = TAU / 8.0;
        for n in 0..4 {
            let k = q.phase[n] / quantum;
            assert_close(k, k.round(), 1e-9);
        }
        // True time delay: only the offsets are phases.
        let mut t = b.clone().with_phase_bits(Some(3)).unwrap();
        t.steer(Direction::from_theta_phi(0.3, 0.0), SteeringMode::TrueTimeDelay);
        assert_close(t.phase[0], 0.0, 1e-12);
        assert_close(t.phase[2], quantum * 0.0 + (0.3f64 / quantum).round() * quantum, 1e-12);
        assert_eq!(t.positions().len(), 4);
        assert_close(t.positions()[1].x - t.positions()[0].x, lambda() / 2.0, 1e-12);
        assert_eq!(t.failed(), &[false; 4]);
    }

    #[test]
    fn invalid_excitations_are_rejected() {
        let a = ula(2);
        assert!(a.clone().with_excitation(&[1.0, 1.0], &[0.0]).is_err());
        assert!(a.clone().with_excitation(&[1.0, 1.0], &[0.0, f64::NAN]).is_err());
        assert!(a.clone().with_excitation(&[1.0], &[0.0, 0.0]).is_err());
        assert!(a.with_excitation(&[0.0, 0.0], &[0.0, 0.0]).is_err());
    }
}
