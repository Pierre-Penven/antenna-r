//! Satellite antenna patterns: ITU-R S.672-4 (GSO, single feed) and
//! S.1528-1 (non-GSO, `recommends 1.2`, 1.3 and 1.4).

use std::f64::consts::PI;

use super::{BeamEllipse, SideLobeEnvelope, off_axis_deg, plane_rad};
use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite, require_level_db, require_positive};
use crate::math::bessel_jn;
use crate::pattern::{DEFAULT_FLOOR_DBI, GainPattern};

/// Largest half beamwidth accepted by the satellite envelopes (degrees).
const MAX_HALF_BEAMWIDTH_DEG: f64 = 90.0;
/// Main-lobe exponent `α` of S.672-4 Table 1.
const S672_ALPHA: f64 = 2.0;
/// Main-lobe exponent `α` of S.1528-1 Table 1.
const S1528_ALPHA: f64 = 1.5;

fn half_beamwidth_ellipse(half_major_rad: f64, half_minor_rad: f64, orientation_rad: f64) -> Result<BeamEllipse> {
    BeamEllipse::new(
        ("half_major_rad", "half_minor_rad"),
        half_major_rad,
        half_minor_rad,
        orientation_rad,
        MAX_HALF_BEAMWIDTH_DEG,
    )
}

/// Near-in side-lobe level `LN` of ITU-R S.672-4 Table 1. `LN = −30 dB` is
/// not offered: the recommendation leaves `a` and `α` for it under study.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum S672SideLobe {
    /// `LN = −20 dB`, `a = 2.58 √(1 − log z)`.
    Minus20,
    /// `LN = −25 dB`, `a = 2.58 √(1 − 0.8 log z)`.
    Minus25,
}

impl S672SideLobe {
    fn level_db(self) -> f64 {
        match self {
            Self::Minus20 => -20.0,
            Self::Minus25 => -25.0,
        }
    }

    fn log_z_coefficient(self) -> f64 {
        match self {
            Self::Minus20 => 1.0,
            Self::Minus25 => 0.8,
        }
    }
}

/// `a = 2.58 √(1 − k log z)` of S.672-4 and S.1528-1 Table 1.
fn table1_a(k: f64, log_z: f64) -> f64 {
    2.58 * (1.0 - k * log_z).sqrt()
}

/// ITU-R S.672-4 `recommends 1`: design-objective pattern of single-feed
/// circular or elliptical beam GSO satellite antennas (ψ in degrees):
///
/// ```text
/// G = Gm − 3 (ψ/ψb)²              ψ ≤ a ψb            (1)
///   = Gm + LN + 20 log z          a ψb < ψ ≤ 0.5 b ψb  (2a)
///   = Gm + LN                     0.5 b ψb < ψ ≤ b ψb  (2b)
///   = X − 25 log ψ                b ψb < ψ ≤ Y         (3)
///   = LF = 0 dBi                  Y < ψ ≤ 90°          (4a)
///   = LB                          90° < ψ ≤ 180°       (4b)
/// X = Gm + LN + 25 log(b ψb),  Y = b ψb 10^(0.04 (Gm + LN − LF))
/// LB = max(15 + LN + 0.25 Gm + 5 log z, 0) dBi
/// a = 2.58 √(1 − log z) (LN = −20 dB) or 2.58 √(1 − 0.8 log z) (−25 dB),
/// b = 6.32
/// ```
///
/// `ψb` is half the −3 dB beamwidth in the plane of interest: an ellipse
/// with semi-axes `ψb,major` (along azimuth `orientation`) and `ψb,minor`;
/// `z = ψb,major / ψb,minor`. Equation (1) is also used inside the coverage
/// (`ψ < ψb`), where the recommendation sets no objective.
/// `recommends 2` (multiple-feed shaped beams) needs the coverage contour
/// and is not modelled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuS672 {
    envelope: SideLobeEnvelope,
}

impl ItuS672 {
    /// Circular beam with half −3 dB beamwidth `half_beamwidth_rad`.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::{ItuS672, S672SideLobe}};
    /// let p = ItuS672::circular(40.0, 1f64.to_radians(), S672SideLobe::Minus20).unwrap();
    /// let g = p.gain_dbi(Direction::from_theta_phi(4f64.to_radians(), 0.0));
    /// assert_eq!(g, 20.0); // Gm + LN on the near side-lobe plateau
    /// ```
    ///
    /// # Errors
    /// As [`ItuS672::elliptical`].
    pub fn circular(peak_gain_dbi: f64, half_beamwidth_rad: f64, side_lobe: S672SideLobe) -> Result<Self> {
        Self::elliptical(peak_gain_dbi, half_beamwidth_rad, half_beamwidth_rad, 0.0, side_lobe)
    }

    /// Elliptical beam: half −3 dB beamwidths along and across the major
    /// axis, which lies at azimuth `orientation_rad`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain or orientation, half-beamwidths outside
    /// `(0, 90°)`, a minor axis larger than the major one, or an axis ratio
    /// for which Table 1 gives no positive `a`.
    pub fn elliptical(
        peak_gain_dbi: f64,
        half_major_rad: f64,
        half_minor_rad: f64,
        orientation_rad: f64,
        side_lobe: S672SideLobe,
    ) -> Result<Self> {
        let ellipse = half_beamwidth_ellipse(half_major_rad, half_minor_rad, orientation_rad)?;
        let k = side_lobe.log_z_coefficient();
        Ok(Self {
            envelope: SideLobeEnvelope::new(peak_gain_dbi, ellipse, side_lobe.level_db(), S672_ALPHA, |log_z| {
                table1_a(k, log_z)
            })?,
        })
    }

    /// Half −3 dB beamwidth `ψb` (degrees) in the plane `φ = phi_rad`.
    #[must_use]
    pub fn half_beamwidth_deg(&self, phi_rad: f64) -> f64 {
        self.envelope.half_beamwidth.width_deg(phi_rad)
    }

    /// Gain at off-axis angle `psi_deg` in a plane where the half
    /// beamwidth is `psi_b_deg`.
    #[must_use]
    pub fn gain_at_deg(&self, psi_deg: f64, psi_b_deg: f64) -> f64 {
        self.envelope.gain_at_deg(psi_deg, psi_b_deg)
    }
}

impl GainPattern for ItuS672 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.envelope.gain_dbi(dir)
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.envelope.peak_gain_dbi
    }
}

/// Near-in side-lobe level `LN` of ITU-R S.1528-1 Table 1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum S1528SideLobe {
    /// `LN = −15 dB`, `a = 2.58 √(1 − 1.4 log z)`.
    Minus15,
    /// `LN = −20 dB`, `a = 2.58 √(1 − 1.0 log z)`.
    Minus20,
    /// `LN = −25 dB`, `a = 2.58 √(1 − 0.6 log z)`.
    Minus25,
    /// `LN = −30 dB`, `a = 2.58 √(1 − 0.4 log z)` (values still under
    /// study in the recommendation).
    Minus30,
}

impl S1528SideLobe {
    fn level_db(self) -> f64 {
        match self {
            Self::Minus15 => -15.0,
            Self::Minus20 => -20.0,
            Self::Minus25 => -25.0,
            Self::Minus30 => -30.0,
        }
    }

    fn log_z_coefficient(self) -> f64 {
        match self {
            Self::Minus15 => 1.4,
            Self::Minus20 => 1.0,
            Self::Minus25 => 0.6,
            Self::Minus30 => 0.4,
        }
    }
}

/// ITU-R S.1528-1 `recommends 1.2`: pattern of multiple-beam non-GSO FSS
/// satellite antennas below 30 GHz, for design objectives or interference
/// analysis. Same six segments as [`ItuS672`], with Table 1 of S.1528-1:
/// `α = 1.5`, `a = 2.58 √(1 − k log z)` with `k` = 1.4, 1.0, 0.6, 0.4 for
/// `LN` = −15, −20, −25, −30 dB, and `b = 6.32`:
///
/// ```text
/// G = Gm − 3 (ψ/ψb)^1.5           0 < ψ ≤ a ψb
///   = Gm + LN + 20 log z          a ψb < ψ ≤ 0.5 b ψb
///   = Gm + LN                     0.5 b ψb < ψ ≤ b ψb
///   = X − 25 log ψ                b ψb < ψ ≤ Y
///   = LF = 0 dBi                  Y < ψ ≤ 90°
///   = LB                          90° < ψ ≤ 180°
/// ```
///
/// When the beamwidth is unknown the recommendation suggests
/// `ψb = √1200 / (D/λ)` degrees for the minor axis (see
/// [`ItuS1528::half_beamwidth_from_d_over_lambda_rad`]); the major axis is
/// `z` times larger. Measured patterns should be preferred when available
/// (`recommends 1.1`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuS1528 {
    envelope: SideLobeEnvelope,
}

impl ItuS1528 {
    /// Circular beam with half −3 dB beamwidth `half_beamwidth_rad`.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::{ItuS1528, S1528SideLobe}};
    /// let p = ItuS1528::circular(40.0, 1f64.to_radians(), S1528SideLobe::Minus20).unwrap();
    /// let g = p.gain_dbi(Direction::from_theta_phi(1f64.to_radians(), 0.0));
    /// assert!((g - 37.0).abs() < 1e-12); // Gm − 3 at ψb
    /// ```
    ///
    /// # Errors
    /// As [`ItuS1528::elliptical`].
    pub fn circular(peak_gain_dbi: f64, half_beamwidth_rad: f64, side_lobe: S1528SideLobe) -> Result<Self> {
        Self::elliptical(peak_gain_dbi, half_beamwidth_rad, half_beamwidth_rad, 0.0, side_lobe)
    }

    /// Elliptical beam: half −3 dB beamwidths along and across the major
    /// axis, which lies at azimuth `orientation_rad`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain or orientation, half-beamwidths outside
    /// `(0, 90°)`, a minor axis larger than the major one, or an axis ratio
    /// for which Table 1 gives no positive `a`.
    pub fn elliptical(
        peak_gain_dbi: f64,
        half_major_rad: f64,
        half_minor_rad: f64,
        orientation_rad: f64,
        side_lobe: S1528SideLobe,
    ) -> Result<Self> {
        let ellipse = half_beamwidth_ellipse(half_major_rad, half_minor_rad, orientation_rad)?;
        let k = side_lobe.log_z_coefficient();
        Ok(Self {
            envelope: SideLobeEnvelope::new(peak_gain_dbi, ellipse, side_lobe.level_db(), S1528_ALPHA, |log_z| {
                table1_a(k, log_z)
            })?,
        })
    }

    /// Minor-axis half beamwidth `√1200 / (D/λ)` degrees, returned in
    /// radians, for use when the actual value is unknown.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-positive `D/λ`.
    pub fn half_beamwidth_from_d_over_lambda_rad(d_over_lambda: f64) -> Result<f64> {
        let dl = require_positive("d_over_lambda", d_over_lambda)?;
        Ok((1200f64.sqrt() / dl).to_radians())
    }

    /// Half −3 dB beamwidth `ψb` (degrees) in the plane `φ = phi_rad`.
    #[must_use]
    pub fn half_beamwidth_deg(&self, phi_rad: f64) -> f64 {
        self.envelope.half_beamwidth.width_deg(phi_rad)
    }

    /// Gain at off-axis angle `psi_deg` in a plane where the half
    /// beamwidth is `psi_b_deg`.
    #[must_use]
    pub fn gain_at_deg(&self, psi_deg: f64, psi_b_deg: f64) -> f64 {
        self.envelope.gain_at_deg(psi_deg, psi_b_deg)
    }
}

impl GainPattern for ItuS1528 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.envelope.gain_dbi(dir)
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.envelope.peak_gain_dbi
    }
}

/// Orbit class of ITU-R S.1528-1 `recommends 1.3`, which fixes the
/// main-beam / near-in side-lobe cross point `Ls`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum S1528Orbit {
    /// Medium Earth orbit: `Ls = −12 dB`, `Y = 2 ψb`.
    Meo,
    /// Low Earth orbit: `Ls = −6.75 dB`, `Y = 1.5 ψb`.
    Leo,
}

impl S1528Orbit {
    fn cross_point_db(self) -> f64 {
        match self {
            Self::Meo => -12.0,
            Self::Leo => -6.75,
        }
    }
}

/// ITU-R S.1528-1 `recommends 1.3`: non-GSO satellite pattern for
/// `D/λ < 35` (ψ in degrees):
///
/// ```text
/// G = Gm − 3 (ψ/ψb)²            ψ < Y
///   = Gm + Ls − 25 log(ψ/Y)     Y ≤ ψ < Z
///   = LF                        Z ≤ ψ ≤ 180°
/// Y = ψb (|Ls|/3)^½,  Z = Y · 10^(0.04 (Gm + Ls − LF))
/// ```
///
/// `ψb` is half the −3 dB beamwidth in the plane of interest at the largest
/// off-axis (scan) angle, `LF` the far side-lobe level (about 0 dBi for
/// ideal patterns). The recommendation also writes the side lobes as
/// `20 log(D/λ) + 3.5 − 25 log(ψ/ψb)` (MEO) and
/// `20 log(D/λ) + 5.62 − 25 log(ψ/ψb)` (LEO); both equal the form above when
/// `Gm ≈ 20 log(D/λ) + 7.97`, so the crate uses the `Gm` form.
/// Equation (1) is also used for `ψ < ψb`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuS1528Ngso {
    peak_gain_dbi: f64,
    half_beamwidth: BeamEllipse,
    cross_point_db: f64,
    y_over_psi_b: f64,
    far_dbi: f64,
}

impl ItuS1528Ngso {
    /// Circular beam.
    ///
    /// ```
    /// use antenna_r::parabolic::itu::{ItuS1528Ngso, S1528Orbit};
    /// // Annex 1 MEO example: ψb = 1.6°, Gm = 35 dBi, LF = 3 dBi.
    /// let p = ItuS1528Ngso::circular(35.0, 1.6f64.to_radians(), S1528Orbit::Meo, 3.0).unwrap();
    /// assert!((p.gain_at_deg(10.0, 1.6) - (35.6 - 25.0)).abs() < 0.05);
    /// ```
    ///
    /// # Errors
    /// As [`ItuS1528Ngso::elliptical`].
    pub fn circular(
        peak_gain_dbi: f64,
        half_beamwidth_rad: f64,
        orbit: S1528Orbit,
        far_level_dbi: f64,
    ) -> Result<Self> {
        Self::elliptical(
            peak_gain_dbi,
            half_beamwidth_rad,
            half_beamwidth_rad,
            0.0,
            orbit,
            far_level_dbi,
        )
    }

    /// Elliptical beam: half −3 dB beamwidths along and across the major
    /// axis, which lies at azimuth `orientation_rad`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain, orientation or far level, half-beamwidths
    /// outside `(0, 90°)` or a minor axis larger than the major one, or a far
    /// level not below `Gm + Ls`.
    pub fn elliptical(
        peak_gain_dbi: f64,
        half_major_rad: f64,
        half_minor_rad: f64,
        orientation_rad: f64,
        orbit: S1528Orbit,
        far_level_dbi: f64,
    ) -> Result<Self> {
        let gm = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let half_beamwidth = half_beamwidth_ellipse(half_major_rad, half_minor_rad, orientation_rad)?;
        let lf = require_finite("far_level_dbi", far_level_dbi)?;
        let ls = orbit.cross_point_db();
        if lf >= gm + ls {
            return Err(invalid(
                "far_level_dbi",
                format!("must be below Gm + Ls = {} dBi, got {lf}", gm + ls),
            ));
        }
        Ok(Self {
            peak_gain_dbi: gm,
            half_beamwidth,
            cross_point_db: ls,
            y_over_psi_b: (-ls / 3.0).sqrt(),
            far_dbi: lf,
        })
    }

    /// Half −3 dB beamwidth `ψb` (degrees) in the plane `φ = phi_rad`.
    #[must_use]
    pub fn half_beamwidth_deg(&self, phi_rad: f64) -> f64 {
        self.half_beamwidth.width_deg(phi_rad)
    }

    /// Gain at off-axis angle `psi_deg` in a plane where the half
    /// beamwidth is `psi_b_deg`.
    #[must_use]
    pub fn gain_at_deg(&self, psi_deg: f64, psi_b_deg: f64) -> f64 {
        let y = psi_b_deg * self.y_over_psi_b;
        let gm_ls = self.peak_gain_dbi + self.cross_point_db;
        if psi_deg < y {
            let r = psi_deg / psi_b_deg;
            self.peak_gain_dbi - 3.0 * r * r
        } else if psi_deg < y * 10f64.powf(0.04 * (gm_ls - self.far_dbi)) {
            gm_ls - 25.0 * (psi_deg / y).log10()
        } else if psi_deg.is_nan() {
            f64::NAN
        } else {
            self.far_dbi
        }
    }
}

impl GainPattern for ItuS1528Ngso {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let psi_b = self.half_beamwidth.width_deg(plane_rad(dir));
        self.gain_at_deg(off_axis_deg(dir), psi_b)
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

/// Largest number of side lobes `l` accepted by [`ItuS1528Taylor`].
pub const S1528_TAYLOR_MAX_LOBES: u32 = 64;
/// `|1 − (u/(π μi))²|` below which the removable singularity of the
/// S.1528 `recommends 1.4` function is replaced by its limit.
const TAYLOR_SINGULAR_EPS: f64 = 1e-10;
/// Newton iterations for the zeros of `J1` (quadratic convergence from
/// the `McMahon` estimate, error below 1e-3).
const J1_ZERO_NEWTON_STEPS: usize = 6;

/// `k`-th positive zero of `J1` (`k ≥ 1`), by Newton's method from
/// the `McMahon` asymptotic estimate.
fn bessel_j1_zero(k: u32) -> f64 {
    let beta = (f64::from(k) + 0.25) * PI;
    let mut x = beta - 3.0 / (8.0 * beta);
    for _ in 0..J1_ZERO_NEWTON_STEPS {
        let j1 = bessel_jn(1, x);
        x -= j1 / (bessel_jn(0, x) - j1 / x);
    }
    x
}

/// One factor of the S.1528 `recommends 1.4` product.
#[derive(Debug, Clone, Copy, PartialEq)]
struct TaylorFactor {
    /// `π μi`, the `i`-th zero of `J1`.
    zero: f64,
    /// `1 / (π² σ² (A² + (i − ½)²))`.
    numerator_coeff: f64,
    /// `lim J1(u) / (1 − (u/(π μi))²)` at `u = π μi`, i.e. `−π μi J0(π μi) / 2`.
    singular_limit: f64,
}

/// ITU-R S.1528-1 `recommends 1.4`: analytical non-GSO satellite pattern
/// based on a circular Taylor illumination, with symmetric side-lobe peaks:
///
/// ```text
/// G(u) = Gmax + 20 log | 2 J1(u)/u · Π_{i=1}^{l−1} [1 − u² / (π² σ² (A² + (i − ½)²))]
///                                            / [1 − (u / (π μi))²] |
/// A = arccosh(10^(SLR/20)) / π,  σ = μl / √(A² + (l − ½)²)
/// u = (π/λ) √((Lr sinθ cosφ)² + (Lt sinθ sinφ)²)
/// ```
///
/// `μi` are the zeros of `J1(π μ)`, `SLR` the first side-lobe ratio (dB),
/// `l` the number of side lobes modelled, `Lr` and `Lt` the radial and
/// transverse sizes of the effective radiating area (along antenna X and Y).
/// The rendered text prints `μ1` in `σ`; the crate uses `μl`, which matches
/// the Annex 2 example (`SLR = 20 dB`, `l = 4`: `A = 0.95277`,
/// `σ = 1.1692`). The removable `0/0` at `u = π μi` uses its limit. Nulls
/// are clamped to `floor_dbi` (the recommendation leaves their depth to
/// further study); the back hemisphere is at `floor_dbi`.
#[derive(Debug, Clone, PartialEq)]
pub struct ItuS1528Taylor {
    peak_gain_dbi: f64,
    floor_dbi: f64,
    k_radial: f64,
    k_transverse: f64,
    a: f64,
    sigma: f64,
    factors: Box<[TaylorFactor]>,
}

impl ItuS1528Taylor {
    /// Pattern from the peak gain, side-lobe ratio `slr_db`, number of
    /// side lobes `lobes` (`l`), radiating sizes `lr_m`, `lt_m` and the
    /// frequency.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::ItuS1528Taylor};
    /// let p = ItuS1528Taylor::new(40.0, 20.0, 4, 1.0, 1.0, 12e9).unwrap();
    /// assert!((p.taylor_a() - 0.95277).abs() < 1e-5);
    /// assert!((p.taylor_sigma() - 1.1692).abs() < 1e-4);
    /// assert_eq!(p.gain_dbi(Direction::BORESIGHT), 40.0);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain, a non-positive `SLR`, sizes or frequency, or
    /// `lobes` outside `1..=S1528_TAYLOR_MAX_LOBES`.
    pub fn new(peak_gain_dbi: f64, slr_db: f64, lobes: u32, lr_m: f64, lt_m: f64, freq_hz: f64) -> Result<Self> {
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let slr = require_positive("slr_db", slr_db)?;
        if !(1..=S1528_TAYLOR_MAX_LOBES).contains(&lobes) {
            return Err(invalid(
                "lobes",
                format!("must be in 1..={S1528_TAYLOR_MAX_LOBES}, got {lobes}"),
            ));
        }
        let lr = require_positive("lr_m", lr_m)?;
        let lt = require_positive("lt_m", lt_m)?;
        let lambda = rf_r::fspl::wavelength_m(require_positive("freq_hz", freq_hz)?);
        let a = rf_r::db::amplitude_from_db(slr).acosh() / PI;
        let a2 = a * a;
        let half = |i: u32| f64::from(i) - 0.5;
        let sigma = bessel_j1_zero(lobes) / PI / half(lobes).mul_add(half(lobes), a2).sqrt();
        let factors = (1..lobes)
            .map(|i| {
                let zero = bessel_j1_zero(i);
                TaylorFactor {
                    zero,
                    numerator_coeff: 1.0 / ((PI * sigma).powi(2) * half(i).mul_add(half(i), a2)),
                    singular_limit: -zero * bessel_jn(0, zero) / 2.0,
                }
            })
            .collect();
        Ok(Self {
            peak_gain_dbi: gmax,
            floor_dbi: DEFAULT_FLOOR_DBI,
            k_radial: PI * lr / lambda,
            k_transverse: PI * lt / lambda,
            a,
            sigma,
            factors,
        })
    }

    /// Replaces the floor (default [`DEFAULT_FLOOR_DBI`]).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for `NaN` or `+∞`.
    pub fn with_floor_dbi(mut self, floor_dbi: f64) -> Result<Self> {
        self.floor_dbi = require_level_db("floor_dbi", floor_dbi)?;
        Ok(self)
    }

    /// Parameter `A` of the illumination.
    #[must_use]
    pub fn taylor_a(&self) -> f64 {
        self.a
    }

    /// Parameter `σ` of the illumination.
    #[must_use]
    pub fn taylor_sigma(&self) -> f64 {
        self.sigma
    }

    /// Normalised field `E(u)` (1 at `u = 0`).
    #[must_use]
    pub fn field(&self, u: f64) -> f64 {
        if u == 0.0 {
            return 1.0;
        }
        let mut product = 1.0;
        let mut singular = None;
        for f in &self.factors {
            product *= (u * u).mul_add(-f.numerator_coeff, 1.0);
            let r = u / f.zero;
            let den = r.mul_add(-r, 1.0);
            if den.abs() < TAYLOR_SINGULAR_EPS {
                singular = Some(f.singular_limit);
            } else {
                product /= den;
            }
        }
        2.0 * singular.unwrap_or_else(|| bessel_jn(1, u)) / u * product
    }

    /// Argument `u` towards `dir` (front hemisphere).
    #[must_use]
    pub fn u_parameter(&self, dir: Direction) -> f64 {
        (self.k_radial * dir.u()).hypot(self.k_transverse * dir.v())
    }
}

impl GainPattern for ItuS1528Taylor {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let w = dir.w();
        if w.is_nan() {
            return f64::NAN;
        }
        if w < 0.0 {
            return self.floor_dbi;
        }
        let field = self.field(self.u_parameter(dir));
        (self.peak_gain_dbi + rf_r::db::amplitude_to_db(field.abs())).max(self.floor_dbi)
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;
    use std::f64::consts::FRAC_PI_2;

    fn at_deg(p: &impl GainPattern, deg: f64) -> f64 {
        p.gain_dbi(Direction::from_theta_phi(deg.to_radians(), 0.3))
    }

    #[test]
    fn s672_circular_segments() {
        let gm = 40.0;
        let p = ItuS672::circular(gm, 1f64.to_radians(), S672SideLobe::Minus20).unwrap();
        assert_eq!(at_deg(&p, 0.0), gm);
        assert_close(at_deg(&p, 1.0), gm - 3.0, 1e-9);
        assert_close(at_deg(&p, 2.58), gm - 3.0 * 2.58 * 2.58, 1e-9);
        // z = 1: (2a) and (2b) both give Gm + LN.
        assert_close(at_deg(&p, 3.0), gm - 20.0, 1e-12);
        assert_close(at_deg(&p, 6.0), gm - 20.0, 1e-12);
        let x = gm - 20.0 + 25.0 * 6.32f64.log10();
        assert_close(at_deg(&p, 10.0), x - 25.0, 1e-9);
        let y = 6.32 * 10f64.powf(0.04 * (gm - 20.0));
        assert_close(at_deg(&p, y * 0.999), x - 25.0 * (y * 0.999).log10(), 1e-9);
        assert_eq!(at_deg(&p, 60.0), 0.0);
        assert_eq!(at_deg(&p, 90.0), 0.0);
        // LB = max(15 − 20 + 10 + 0, 0) = 5 dBi.
        assert_close(at_deg(&p, 120.0), 5.0, 1e-12);
        assert_eq!(p.peak_gain_dbi(), gm);
        assert!(p.gain_at_deg(f64::NAN, 1.0).is_nan());
    }

    #[test]
    fn s672_minus25_and_low_gain_back_lobe() {
        let p = ItuS672::circular(30.0, 2f64.to_radians(), S672SideLobe::Minus25).unwrap();
        assert_close(at_deg(&p, 8.0), 5.0, 1e-12);
        // LB = max(15 − 25 + 7.5, 0) = 0
        assert_eq!(at_deg(&p, 150.0), 0.0);
    }

    #[test]
    fn s672_segment_3_reaches_lf_at_y() {
        let gm = 45.0;
        let p = ItuS672::circular(gm, 0.5f64.to_radians(), S672SideLobe::Minus20).unwrap();
        let y = 6.32 * 0.5 * 10f64.powf(0.04 * (gm - 20.0));
        assert_close(p.gain_at_deg(y, 0.5), 0.0, 1e-9);
    }

    #[test]
    fn s672_elliptical_beam_uses_plane_of_interest() {
        let p = ItuS672::elliptical(35.0, 2f64.to_radians(), 1f64.to_radians(), 0.0, S672SideLobe::Minus20).unwrap();
        assert_close(p.half_beamwidth_deg(0.0), 2.0, 1e-12);
        assert_close(p.half_beamwidth_deg(FRAC_PI_2), 1.0, 1e-12);
        // Along the major axis ψ = ψb gives Gm − 3.
        let g = p.gain_dbi(Direction::from_theta_phi(2f64.to_radians(), 0.0));
        assert_close(g, 32.0, 1e-9);
        // z = 2, a = 2.58 √(1 − log 2) ≈ 2.157 (Table 1 of S.672-4).
        let a = 2.58 * (1.0 - 2f64.log10()).sqrt();
        let main = p.gain_dbi(Direction::from_theta_phi(2.1f64.to_radians(), FRAC_PI_2));
        assert_close(main, 35.0 - 3.0 * 2.1 * 2.1, 1e-9);
        assert!(2.1 < a && a < 2.2);
        // On the minor axis (ψb = 1°), 2.8° lies on (2a): Gm + LN + 20 log z.
        let g2a = p.gain_dbi(Direction::from_theta_phi(2.8f64.to_radians(), FRAC_PI_2));
        assert_close(g2a, 35.0 - 20.0 + 20.0 * 2f64.log10(), 1e-9);
    }

    #[test]
    fn s672_rejects_invalid_parameters() {
        let ok = 1f64.to_radians();
        assert!(matches!(
            ItuS672::circular(f64::NAN, ok, S672SideLobe::Minus20),
            Err(AntennaError::InvalidParameter {
                name: "peak_gain_dbi",
                ..
            })
        ));
        assert!(matches!(
            ItuS672::circular(40.0, 0.0, S672SideLobe::Minus20),
            Err(AntennaError::InvalidParameter {
                name: "half_major_rad",
                ..
            })
        ));
        assert!(matches!(
            ItuS672::circular(40.0, 2.0, S672SideLobe::Minus20),
            Err(AntennaError::InvalidParameter {
                name: "half_major_rad",
                ..
            })
        ));
        assert!(matches!(
            ItuS672::elliptical(40.0, ok, 2.0 * ok, 0.0, S672SideLobe::Minus20),
            Err(AntennaError::InvalidParameter {
                name: "half_minor_rad",
                ..
            })
        ));
        // z = 20: 1 − log z < 0, no real a.
        assert!(matches!(
            ItuS672::elliptical(40.0, 20.0 * ok, ok, 0.0, S672SideLobe::Minus20),
            Err(AntennaError::InvalidParameter {
                name: "half_minor_rad",
                ..
            })
        ));
        assert!(matches!(
            ItuS672::elliptical(40.0, ok, ok, f64::NAN, S672SideLobe::Minus25),
            Err(AntennaError::InvalidParameter {
                name: "orientation_rad",
                ..
            })
        ));
    }

    #[test]
    fn s1528_envelope_uses_alpha_1_5_and_table_1() {
        let gm = 40.0;
        let p = ItuS1528::circular(gm, 1f64.to_radians(), S1528SideLobe::Minus20).unwrap();
        assert_eq!(at_deg(&p, 0.0), gm);
        assert_close(at_deg(&p, 2.0), gm - 3.0 * 2f64.powf(1.5), 1e-9);
        assert_close(at_deg(&p, 2.58), gm - 3.0 * 2.58f64.powf(1.5), 1e-9);
        assert_close(at_deg(&p, 3.0), gm - 20.0, 1e-12);
        assert_close(at_deg(&p, 6.32), gm - 20.0, 1e-12);
        assert_close(at_deg(&p, 10.0), gm - 20.0 + 25.0 * 6.32f64.log10() - 25.0, 1e-9);
        assert_eq!(at_deg(&p, 80.0), 0.0);
        assert_close(at_deg(&p, 170.0), 5.0, 1e-12);
        assert_eq!(p.peak_gain_dbi(), gm);
        assert_close(p.half_beamwidth_deg(0.0), 1.0, 1e-12);
        assert!(p.gain_at_deg(f64::NAN, 1.0).is_nan());
        for (sl, ln) in [
            (S1528SideLobe::Minus15, -15.0),
            (S1528SideLobe::Minus25, -25.0),
            (S1528SideLobe::Minus30, -30.0),
        ] {
            let q = ItuS1528::circular(gm, 1f64.to_radians(), sl).unwrap();
            assert_close(q.gain_at_deg(4.0, 1.0), gm + ln, 1e-12);
        }
    }

    #[test]
    fn s1528_elliptical_a_depends_on_side_lobe_level() {
        let deg = 1f64.to_radians();
        // z = 2: a = 2.58 √(1 − 1.4 log 2) ≈ 1.96 for LN = −15 dB.
        let p = ItuS1528::elliptical(35.0, 2.0 * deg, deg, 0.0, S1528SideLobe::Minus15).unwrap();
        let a = 2.58 * (1.0 - 1.4 * 2f64.log10()).sqrt();
        assert_close(p.gain_at_deg(0.99 * a, 1.0), 35.0 - 3.0 * (0.99 * a).powf(1.5), 1e-9);
        assert_close(p.gain_at_deg(1.01 * a, 1.0), 35.0 - 15.0 + 20.0 * 2f64.log10(), 1e-9);
        // z = 6: 1 − 1.4 log 6 < 0.
        assert!(ItuS1528::elliptical(35.0, 6.0 * deg, deg, 0.0, S1528SideLobe::Minus15).is_err());
        assert!(ItuS1528::elliptical(35.0, 6.0 * deg, deg, 0.0, S1528SideLobe::Minus30).is_ok());
        let psi_b = ItuS1528::half_beamwidth_from_d_over_lambda_rad(20.0).unwrap();
        assert_close(psi_b.to_degrees(), 1200f64.sqrt() / 20.0, 1e-12);
        assert!(ItuS1528::half_beamwidth_from_d_over_lambda_rad(0.0).is_err());
    }

    #[test]
    fn s1528_ngso_matches_annex_1_examples() {
        // MEO: ψb = 1.6°, Gm = 35, Ls = −12, LF = 3 → Y = 3.2°, Z ≈ 20°.
        let meo = ItuS1528Ngso::circular(35.0, 1.6f64.to_radians(), S1528Orbit::Meo, 3.0).unwrap();
        assert_close(meo.gain_at_deg(1.6, 1.6), 32.0, 1e-12);
        assert_close(meo.gain_at_deg(3.2, 1.6), 23.0, 1e-12);
        assert_close(meo.gain_at_deg(10.0, 1.6), 23.0 - 25.0 * (10.0f64 / 3.2).log10(), 1e-12);
        assert_eq!(meo.gain_at_deg(20.3, 1.6), 3.0);
        assert_close(meo.gain_at_deg(20.1, 1.6), 23.0 - 25.0 * (20.1f64 / 3.2).log10(), 1e-12);
        // LEO: Ls = −6.75, LF = 5 → Y = 2.4°, Z ≈ 20.4°, 37.76 − 25 log ψ.
        let leo = ItuS1528Ngso::circular(35.0, 1.6f64.to_radians(), S1528Orbit::Leo, 5.0).unwrap();
        assert_close(leo.gain_at_deg(2.4, 1.6), 28.25, 1e-12);
        assert_close(leo.gain_at_deg(10.0, 1.6), 37.76 - 25.0, 0.01);
        assert_eq!(leo.gain_at_deg(20.5, 1.6), 5.0);
        assert_eq!(at_deg(&leo, 0.0), 35.0);
        assert_eq!(leo.peak_gain_dbi(), 35.0);
        assert!(leo.gain_at_deg(f64::NAN, 1.6).is_nan());
    }

    #[test]
    fn s1528_ngso_elliptical_and_validation() {
        let deg = 1f64.to_radians();
        let p = ItuS1528Ngso::elliptical(35.0, 2.0 * deg, deg, 0.0, S1528Orbit::Meo, 0.0).unwrap();
        assert_close(p.half_beamwidth_deg(FRAC_PI_2), 1.0, 1e-12);
        let g = p.gain_dbi(Direction::from_theta_phi(2f64.to_radians(), 0.0));
        assert_close(g, 32.0, 1e-12);
        assert!(matches!(
            ItuS1528Ngso::circular(35.0, deg, S1528Orbit::Meo, 23.0),
            Err(AntennaError::InvalidParameter {
                name: "far_level_dbi",
                ..
            })
        ));
        assert!(ItuS1528Ngso::circular(f64::NAN, deg, S1528Orbit::Meo, 0.0).is_err());
        assert!(ItuS1528Ngso::circular(35.0, 0.0, S1528Orbit::Meo, 0.0).is_err());
        assert!(ItuS1528Ngso::circular(35.0, deg, S1528Orbit::Leo, f64::NAN).is_err());
    }

    #[test]
    fn j1_zeros_match_the_recommendation() {
        // Annex 2: μ = 1.220, 2.233, 3.238, 4.241.
        for (k, mu) in [(1, 1.2197), (2, 2.2331), (3, 3.2383), (4, 4.2411)] {
            let z = bessel_j1_zero(k);
            assert_close(z / PI, mu, 1e-4);
            assert!(bessel_jn(1, z).abs() < 1e-12);
        }
        assert!(bessel_jn(1, bessel_j1_zero(40)).abs() < 1e-11);
    }

    #[test]
    fn s1528_taylor_with_one_lobe_is_the_uniform_aperture() {
        let p = ItuS1528Taylor::new(40.0, 20.0, 1, 1.0, 1.0, 12e9).unwrap();
        for u in [0.5, 2.0, 5.0] {
            assert_close(p.field(u), 2.0 * bessel_jn(1, u) / u, 1e-15);
        }
        assert_eq!(p.field(0.0), 1.0);
    }

    #[test]
    fn s1528_taylor_side_lobe_ratio_and_symmetry() {
        // Annex 2 example parameters.
        let annex = ItuS1528Taylor::new(40.0, 20.0, 4, 1.0, 1.0, 12e9).unwrap();
        assert_close(annex.taylor_a(), 0.95277, 1e-5);
        assert_close(annex.taylor_sigma(), 1.1692, 1e-4);
        let lambda = rf_r::fspl::wavelength_m(12e9);
        let p = ItuS1528Taylor::new(40.0, 20.0, 4, 10.0 * lambda, 10.0 * lambda, 12e9).unwrap();
        // First side-lobe peak: −20.72 dB at u ≈ 5.05 for SLR = 20 dB and
        // l = 4 (independent SciPy evaluation of the same formula; the
        // Taylor approximation with few lobes lands slightly below SLR).
        let peak_sl = (400..660)
            .map(|i| f64::from(i) * 0.01)
            .map(|u| 20.0 * p.field(u).abs().log10())
            .fold(f64::NEG_INFINITY, f64::max);
        assert_close(peak_sl, -20.72, 0.01);
        let q = ItuS1528Taylor::new(40.0, 30.0, 8, 1.0, 1.0, 12e9).unwrap();
        let q_sl = (450..700)
            .map(|i| f64::from(i) * 0.01)
            .map(|u| 20.0 * q.field(u).abs().log10())
            .fold(f64::NEG_INFINITY, f64::max);
        assert_close(q_sl, -30.24, 0.01);
        // Lr = Lt: rotationally symmetric.
        let a = p.gain_dbi(Direction::from_theta_phi(3f64.to_radians(), 0.0));
        let b = p.gain_dbi(Direction::from_theta_phi(3f64.to_radians(), 1.0));
        assert_close(a, b, 1e-9);
        assert_close(
            p.u_parameter(Direction::from_theta_phi(FRAC_PI_2, 0.0)),
            10.0 * PI,
            1e-9,
        );
    }

    #[test]
    fn s1528_taylor_removable_singularity_is_continuous() {
        let p = ItuS1528Taylor::new(40.0, 25.0, 5, 1.0, 1.0, 12e9).unwrap();
        let z = bessel_j1_zero(2);
        let exact = p.field(z);
        assert_close(exact, p.field(z * (1.0 + 1e-7)), 1e-6);
        assert_close(exact, p.field(z * (1.0 - 1e-7)), 1e-6);
        assert!(exact.abs() > 1e-6);
    }

    #[test]
    fn s1528_taylor_hemispheres_floor_and_nan() {
        let p = ItuS1528Taylor::new(40.0, 20.0, 4, 0.5, 0.3, 12e9)
            .unwrap()
            .with_floor_dbi(-20.0)
            .unwrap();
        assert_eq!(p.gain_dbi(Direction::BORESIGHT), 40.0);
        assert_eq!(p.peak_gain_dbi(), 40.0);
        assert_eq!(p.gain_dbi(Direction::from_theta_phi(2.0, 0.0)), -20.0);
        assert!(p.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
        let g = p.gain_dbi(Direction::from_theta_phi(0.2, 0.0));
        assert!((-20.0..40.0).contains(&g));
        assert!(p.clone().with_floor_dbi(f64::NAN).is_err());
    }

    #[test]
    fn s1528_taylor_validation() {
        let ok = ItuS1528Taylor::new(40.0, 20.0, 4, 1.0, 1.0, 12e9);
        assert!(ok.is_ok());
        for (args, name) in [
            ((f64::NAN, 20.0, 4, 1.0, 1.0, 12e9), "peak_gain_dbi"),
            ((40.0, 0.0, 4, 1.0, 1.0, 12e9), "slr_db"),
            ((40.0, 20.0, 0, 1.0, 1.0, 12e9), "lobes"),
            ((40.0, 20.0, 65, 1.0, 1.0, 12e9), "lobes"),
            ((40.0, 20.0, 4, 0.0, 1.0, 12e9), "lr_m"),
            ((40.0, 20.0, 4, 1.0, -1.0, 12e9), "lt_m"),
            ((40.0, 20.0, 4, 1.0, 1.0, 0.0), "freq_hz"),
        ] {
            let (g, s, l, lr, lt, f) = args;
            assert!(
                matches!(
                    ItuS1528Taylor::new(g, s, l, lr, lt, f),
                    Err(AntennaError::InvalidParameter { name: n, .. }) if n == name
                ),
                "{name}"
            );
        }
    }
}
