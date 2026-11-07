//! Parabolic and circular-aperture antennas: analytic aperture patterns, a
//! Gaussian main beam and (in [`itu`]) ITU-R reference envelopes.
//!
//! Peak gains use `rf_r::antenna::parabolic_gain_dbi`. The aperture models
//! are front-hemisphere models: behind the aperture (`θ ≥ 90°`) they return
//! their floor.

pub mod itu;

use std::f64::consts::PI;

use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite, require_level_db, require_positive};
use crate::math::{ellipse_radius, lambda_n};
use crate::pattern::{DEFAULT_FLOOR_DBI, GainPattern};

/// Validates an efficiency in `(0, 1]`.
pub(crate) fn require_efficiency(name: &'static str, value: f64) -> Result<f64> {
    if value > 0.0 && value <= 1.0 {
        Ok(value)
    } else {
        Err(invalid(name, format!("must be in (0, 1], got {value}")))
    }
}

/// Uniformly illuminated circular aperture:
/// `G(θ) = G0 · [2 J1(x)/x]²` with `x = π D sinθ / λ`.
///
/// `G0 = η (π D / λ)²`. First null at `sinθ = 1.2197 λ/D`, first side lobe
/// −17.57 dB.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformAperture {
    peak_gain_dbi: f64,
    k_radius: f64,
    floor_dbi: f64,
}

impl UniformAperture {
    /// Aperture of `diameter_m` at `freq_hz` with aperture efficiency
    /// `efficiency` in `(0, 1]`.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::UniformAperture};
    /// let dish = UniformAperture::new(1.2, 12e9, 1.0).unwrap();
    /// assert!((dish.peak_gain_dbi() - 43.58).abs() < 0.01);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-positive diameter or frequency, or an efficiency outside `(0, 1]`.
    pub fn new(diameter_m: f64, freq_hz: f64, efficiency: f64) -> Result<Self> {
        let d = require_positive("diameter_m", diameter_m)?;
        let f = require_positive("freq_hz", freq_hz)?;
        let eta = require_efficiency("efficiency", efficiency)?;
        Ok(Self {
            peak_gain_dbi: rf_r::antenna::parabolic_gain_dbi(d, f, eta),
            k_radius: PI * d / rf_r::fspl::wavelength_m(f),
            floor_dbi: DEFAULT_FLOOR_DBI,
        })
    }

    /// Sets the floor (nulls and back hemisphere) in dBi.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// unless finite or `-inf`.
    pub fn with_floor_dbi(mut self, floor_dbi: f64) -> Result<Self> {
        self.floor_dbi = require_level_db("floor_dbi", floor_dbi)?;
        Ok(self)
    }

    /// Normalised aperture parameter `π D / λ`.
    #[must_use]
    pub fn k_radius(&self) -> f64 {
        self.k_radius
    }
}

impl GainPattern for UniformAperture {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        aperture_gain(self.peak_gain_dbi, self.floor_dbi, self.k_radius, dir, |x| {
            lambda_n(0, x)
        })
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

/// Shared evaluation of circular-aperture patterns from a normalised field
/// function of `x = k a sinθ` (`k_radius = π D / λ`).
#[inline]
fn aperture_gain(peak_dbi: f64, floor_dbi: f64, k_radius: f64, dir: Direction, field: impl Fn(f64) -> f64) -> f64 {
    let w = dir.w();
    if w.is_nan() {
        return f64::NAN;
    }
    if w <= 0.0 {
        return floor_dbi;
    }
    let x = k_radius * dir.u().hypot(dir.v());
    (peak_dbi + rf_r::db::amplitude_to_db(field(x).abs())).max(floor_dbi)
}

/// Circular aperture with a parabolic-on-pedestal illumination
/// `f(ρ) = c + (1 − c)(1 − ρ²)^n`, `c = 10^(edge_taper_db / 20)`.
///
/// Normalised field: `E(x) = [c Λ0(x)/2 + (1 − c) Λn(x) / (2(n + 1))] / E(0)`
/// (Silver, *Microwave Antenna Theory and Design*, §6.7). The taper
/// efficiency `η_t = 2 (∫ f ρ dρ)² / ∫ f² ρ dρ` multiplies the other
/// efficiencies in the peak gain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TaperedAperture {
    peak_gain_dbi: f64,
    k_radius: f64,
    n: u32,
    pedestal: f64,
    taper_efficiency: f64,
    floor_dbi: f64,
}

impl TaperedAperture {
    /// Tapered aperture of `diameter_m` at `freq_hz`, taper exponent
    /// `n ∈ {1, 2}`, edge illumination `edge_taper_db ≤ 0` (power, e.g. −10)
    /// and other efficiencies (spill-over, ohmic, …) `other_efficiency`.
    ///
    /// ```
    /// use antenna_r::parabolic::TaperedAperture;
    /// let a = TaperedAperture::new(1.0, 10e9, 1, f64::NEG_INFINITY, 1.0).unwrap();
    /// assert!((a.taper_efficiency() - 0.75).abs() < 1e-12);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for invalid sizes, `n ∉ {1, 2}`, a positive or NaN edge taper, or an
    /// efficiency outside `(0, 1]`.
    pub fn new(diameter_m: f64, freq_hz: f64, n: u32, edge_taper_db: f64, other_efficiency: f64) -> Result<Self> {
        let d = require_positive("diameter_m", diameter_m)?;
        let f = require_positive("freq_hz", freq_hz)?;
        let eta = require_efficiency("other_efficiency", other_efficiency)?;
        if !(n == 1 || n == 2) {
            return Err(invalid("n", format!("taper exponent must be 1 or 2, got {n}")));
        }
        if edge_taper_db.is_nan() || edge_taper_db > 0.0 {
            return Err(invalid("edge_taper_db", format!("must be ≤ 0 dB, got {edge_taper_db}")));
        }
        let c = rf_r::db::amplitude_from_db(edge_taper_db);
        let taper_efficiency = taper_efficiency(n, c);
        Ok(Self {
            peak_gain_dbi: rf_r::antenna::parabolic_gain_dbi(d, f, eta * taper_efficiency),
            k_radius: PI * d / rf_r::fspl::wavelength_m(f),
            n,
            pedestal: c,
            taper_efficiency,
            floor_dbi: DEFAULT_FLOOR_DBI,
        })
    }

    /// Taper (illumination) efficiency in `(0, 1]`.
    #[must_use]
    pub fn taper_efficiency(&self) -> f64 {
        self.taper_efficiency
    }

    /// Sets the floor (nulls and back hemisphere) in dBi.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// unless finite or `-inf`.
    pub fn with_floor_dbi(mut self, floor_dbi: f64) -> Result<Self> {
        self.floor_dbi = require_level_db("floor_dbi", floor_dbi)?;
        Ok(self)
    }

    fn field(&self, x: f64) -> f64 {
        let c = self.pedestal;
        let m = f64::from(self.n + 1);
        let e = 0.5 * c * lambda_n(0, x) + (1.0 - c) * lambda_n(self.n, x) / (2.0 * m);
        e / (0.5 * c + (1.0 - c) / (2.0 * m))
    }
}

/// `η_t = (∫ f ρ dρ)² / (½ ∫ f² ρ dρ)` for `f = c + (1 − c)(1 − ρ²)^n`.
fn taper_efficiency(n: u32, c: f64) -> f64 {
    let nf = f64::from(n);
    let mean = 0.5 * c + (1.0 - c) / (2.0 * (nf + 1.0));
    let power = 0.5 * c * c + c * (1.0 - c) / (nf + 1.0) + (1.0 - c).powi(2) / (2.0 * (2.0 * nf + 1.0));
    mean * mean / (0.5 * power)
}

impl GainPattern for TaperedAperture {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        aperture_gain(self.peak_gain_dbi, self.floor_dbi, self.k_radius, dir, |x| {
            self.field(x)
        })
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

/// Gaussian main beam, circular or elliptical:
/// `G = Gmax − 12 (θ / θ3dB(φ))²`, floored at `floor_dbi`.
///
/// `θ3dB(φ)` is the full −3 dB beamwidth in the plane `φ`: an ellipse with
/// beamwidths `bw_major` along `φ = orientation` and `bw_minor` across it.
/// The model is only meaningful in the main lobe; the floor stands in for
/// the side lobes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianBeam {
    peak_gain_dbi: f64,
    bw_major_rad: f64,
    bw_minor_rad: f64,
    orientation_rad: f64,
    floor_dbi: f64,
}

/// Rule-of-thumb factor of the −3 dB beamwidth `k λ / D` (degrees).
pub const DEFAULT_BEAMWIDTH_FACTOR_DEG: f64 = rf_r::antenna::HPBW_FACTOR_DEG;

impl GaussianBeam {
    /// Circular Gaussian beam with full −3 dB beamwidth `beamwidth_3db_rad`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain or a beamwidth outside `(0, π]`.
    pub fn new(peak_gain_dbi: f64, beamwidth_3db_rad: f64) -> Result<Self> {
        Self::elliptical(peak_gain_dbi, beamwidth_3db_rad, beamwidth_3db_rad, 0.0)
    }

    /// Elliptical Gaussian beam: beamwidths along and across the major axis
    /// at azimuth `orientation_rad`.
    ///
    /// # Errors
    /// As [`GaussianBeam::new`], plus a non-finite orientation.
    pub fn elliptical(peak_gain_dbi: f64, bw_major_rad: f64, bw_minor_rad: f64, orientation_rad: f64) -> Result<Self> {
        let check_bw = |name, bw: f64| {
            if bw > 0.0 && bw <= PI {
                Ok(bw)
            } else {
                Err(invalid(name, format!("must be in (0, π], got {bw}")))
            }
        };
        Ok(Self {
            peak_gain_dbi: require_finite("peak_gain_dbi", peak_gain_dbi)?,
            bw_major_rad: check_bw("bw_major_rad", bw_major_rad)?,
            bw_minor_rad: check_bw("bw_minor_rad", bw_minor_rad)?,
            orientation_rad: require_finite("orientation_rad", orientation_rad)?,
            floor_dbi: DEFAULT_FLOOR_DBI,
        })
    }

    /// Circular beam of a dish: peak from `parabolic_gain_dbi`, beamwidth
    /// `k λ / D` degrees (`k` = [`DEFAULT_BEAMWIDTH_FACTOR_DEG`] typically).
    ///
    /// ```
    /// use antenna_r::{GainPattern, parabolic::GaussianBeam};
    /// let b = GaussianBeam::from_dish(1.2, 12e9, 0.65, 70.0).unwrap();
    /// assert!((b.peak_gain_dbi() - 41.71).abs() < 0.01);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for invalid dish parameters.
    pub fn from_dish(diameter_m: f64, freq_hz: f64, efficiency: f64, factor_deg: f64) -> Result<Self> {
        let d = require_positive("diameter_m", diameter_m)?;
        let f = require_positive("freq_hz", freq_hz)?;
        let eta = require_efficiency("efficiency", efficiency)?;
        let k = require_positive("factor_deg", factor_deg)?;
        let bw = (k * rf_r::fspl::wavelength_m(f) / d).to_radians();
        Self::new(rf_r::antenna::parabolic_gain_dbi(d, f, eta), bw.min(PI))
    }

    /// Sets the floor in dBi.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// unless finite or `-inf`.
    pub fn with_floor_dbi(mut self, floor_dbi: f64) -> Result<Self> {
        self.floor_dbi = require_level_db("floor_dbi", floor_dbi)?;
        Ok(self)
    }

    /// Full −3 dB beamwidth in the plane `φ = phi_rad`.
    #[must_use]
    pub fn beamwidth_rad(&self, phi_rad: f64) -> f64 {
        ellipse_radius(self.bw_major_rad, self.bw_minor_rad, phi_rad - self.orientation_rad)
    }
}

impl GainPattern for GaussianBeam {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let theta = dir.theta_rad();
        let bw = self.beamwidth_rad(dir.v().atan2(dir.u()));
        let r = theta / bw;
        (self.peak_gain_dbi - 12.0 * r * r).max(self.floor_dbi)
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::{beamwidth_rad, cut_features, mean_gain_db_piecewise};
    use crate::error::AntennaError;
    use crate::math::assert_close;
    use std::f64::consts::FRAC_PI_2;

    const D: f64 = 1.0;
    const F: f64 = 30e9; // D/λ ≈ 100

    fn lambda() -> f64 {
        rf_r::fspl::wavelength_m(F)
    }

    #[test]
    fn uniform_peak_matches_rf_r() {
        let a = UniformAperture::new(D, F, 0.7).unwrap();
        assert_close(a.peak_gain_dbi(), rf_r::antenna::parabolic_gain_dbi(D, F, 0.7), 1e-12);
        assert_close(a.gain_dbi(Direction::BORESIGHT), a.peak_gain_dbi(), 1e-12);
        assert_close(a.k_radius(), PI * D / lambda(), 1e-9);
    }

    #[test]
    fn uniform_first_null_and_side_lobe() {
        let a = UniformAperture::new(D, F, 1.0).unwrap().with_floor_dbi(-300.0).unwrap();
        let f = cut_features(&a, 0.0, 0.05, 20001).unwrap();
        assert_close(f.first_null_rad.sin(), 1.219_67 * lambda() / D, 1e-6);
        assert_close(f.peak_sidelobe_rel_db, -17.57, 0.01);
    }

    #[test]
    fn uniform_half_power_beamwidth_is_1_029_lambda_over_d() {
        let a = UniformAperture::new(D, F, 1.0).unwrap();
        let bw = beamwidth_rad(&a, 0.0, 10.0 * 2.0_f64.log10()).unwrap();
        assert_close(bw, 1.028_99 * lambda() / D, 1e-6);
    }

    #[test]
    fn uniform_is_lossless_when_integrated() {
        let a = UniformAperture::new(0.3, 10e9, 1.0).unwrap();
        let mean = mean_gain_db_piecewise(&a, &[0.0, 0.2, FRAC_PI_2], 400, 1);
        // Front-hemisphere aperture: the integral equals η within the
        // small-angle approximation of the aperture model.
        assert_close(mean, 0.0, 0.05);
    }

    #[test]
    fn uniform_back_hemisphere_and_nan() {
        let a = UniformAperture::new(D, F, 1.0).unwrap();
        assert_eq!(a.gain_dbi(Direction::from_theta_phi(2.0, 0.0)), DEFAULT_FLOOR_DBI);
        assert!(a.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
    }

    #[test]
    fn uniform_rejects_invalid_parameters() {
        assert!(matches!(
            UniformAperture::new(0.0, F, 1.0),
            Err(AntennaError::InvalidParameter { name: "diameter_m", .. })
        ));
        assert!(matches!(
            UniformAperture::new(D, -1.0, 1.0),
            Err(AntennaError::InvalidParameter { name: "freq_hz", .. })
        ));
        for eta in [0.0, 1.1, f64::NAN] {
            assert!(matches!(
                UniformAperture::new(D, F, eta),
                Err(AntennaError::InvalidParameter { name: "efficiency", .. })
            ));
        }
        assert!(
            UniformAperture::new(D, F, 1.0)
                .unwrap()
                .with_floor_dbi(f64::NAN)
                .is_err()
        );
    }

    #[test]
    fn taper_efficiency_reference_values() {
        assert_close(taper_efficiency(1, 0.0), 0.75, 1e-12);
        assert_close(taper_efficiency(2, 0.0), 5.0 / 9.0, 1e-12);
        assert_close(taper_efficiency(1, 1.0), 1.0, 1e-12);
        assert_close(taper_efficiency(2, 1.0), 1.0, 1e-12);
    }

    #[test]
    fn tapered_with_full_pedestal_equals_uniform() {
        let t = TaperedAperture::new(D, F, 1, 0.0, 0.8).unwrap();
        let u = UniformAperture::new(D, F, 0.8).unwrap();
        for theta in [0.0, 0.004, 0.01, 0.03] {
            let d = Direction::from_theta_phi(theta, 0.2);
            assert_close(t.gain_dbi(d), u.gain_dbi(d), 1e-9);
        }
    }

    #[test]
    fn tapered_n1_has_lower_side_lobes() {
        // Zero edge illumination, n = 1: first side lobe −24.64 dB.
        let t = TaperedAperture::new(D, F, 1, f64::NEG_INFINITY, 1.0)
            .unwrap()
            .with_floor_dbi(-300.0)
            .unwrap();
        let f = cut_features(&t, 0.0, 0.05, 20001).unwrap();
        assert_close(f.peak_sidelobe_rel_db, -24.64, 0.05);
        assert_close(f.first_null_rad.sin(), 1.6347 * lambda() / D, 2e-5);
        assert_close(t.peak_gain_dbi(), rf_r::antenna::parabolic_gain_dbi(D, F, 0.75), 1e-12);
        assert_close(t.taper_efficiency(), 0.75, 1e-12);
    }

    #[test]
    fn tapered_n2_side_lobe_level() {
        // Zero edge illumination, n = 2: first side lobe −30.61 dB.
        let t = TaperedAperture::new(D, F, 2, f64::NEG_INFINITY, 1.0)
            .unwrap()
            .with_floor_dbi(-300.0)
            .unwrap();
        let f = cut_features(&t, 0.0, 0.05, 20001).unwrap();
        assert_close(f.peak_sidelobe_rel_db, -30.61, 0.05);
    }

    #[test]
    fn tapered_back_hemisphere_and_validation() {
        let t = TaperedAperture::new(D, F, 2, -10.0, 1.0).unwrap();
        assert_eq!(t.gain_dbi(Direction::from_theta_phi(PI, 0.0)), DEFAULT_FLOOR_DBI);
        assert!(t.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
        assert!(matches!(
            TaperedAperture::new(D, F, 3, -10.0, 1.0),
            Err(AntennaError::InvalidParameter { name: "n", .. })
        ));
        assert!(matches!(
            TaperedAperture::new(D, F, 0, -10.0, 1.0),
            Err(AntennaError::InvalidParameter { name: "n", .. })
        ));
        for et in [1.0, f64::NAN] {
            assert!(matches!(
                TaperedAperture::new(D, F, 1, et, 1.0),
                Err(AntennaError::InvalidParameter {
                    name: "edge_taper_db",
                    ..
                })
            ));
        }
        assert!(TaperedAperture::new(-1.0, F, 1, -10.0, 1.0).is_err());
        assert!(TaperedAperture::new(D, 0.0, 1, -10.0, 1.0).is_err());
        assert!(TaperedAperture::new(D, F, 1, -10.0, 0.0).is_err());
        assert!(t.with_floor_dbi(f64::INFINITY).is_err());
    }

    #[test]
    fn gaussian_follows_its_formula() {
        let g = GaussianBeam::new(40.0, 0.02).unwrap();
        assert_eq!(g.gain_dbi(Direction::BORESIGHT), 40.0);
        assert_close(g.gain_dbi(Direction::from_theta_phi(0.01, 1.0)), 37.0, 1e-12);
        assert_eq!(g.gain_dbi(Direction::from_theta_phi(1.0, 0.0)), DEFAULT_FLOOR_DBI);
        let floored = g.with_floor_dbi(0.0).unwrap();
        assert_eq!(floored.gain_dbi(Direction::from_theta_phi(1.0, 0.0)), 0.0);
        assert_eq!(floored.peak_gain_dbi(), 40.0);
    }

    #[test]
    fn gaussian_elliptical_beamwidths() {
        let g = GaussianBeam::elliptical(30.0, 0.04, 0.02, FRAC_PI_2).unwrap();
        assert_close(g.beamwidth_rad(FRAC_PI_2), 0.04, 1e-15);
        assert_close(g.beamwidth_rad(0.0), 0.02, 1e-15);
        assert_close(g.gain_dbi(Direction::from_theta_phi(0.02, FRAC_PI_2)), 27.0, 1e-12);
        assert_close(g.gain_dbi(Direction::from_theta_phi(0.01, 0.0)), 27.0, 1e-12);
    }

    #[test]
    fn gaussian_from_dish_uses_rule_of_thumb() {
        let g = GaussianBeam::from_dish(D, F, 0.6, DEFAULT_BEAMWIDTH_FACTOR_DEG).unwrap();
        assert_close(g.beamwidth_rad(0.0), (70.0 * lambda() / D).to_radians(), 1e-12);
        assert_close(g.peak_gain_dbi(), rf_r::antenna::parabolic_gain_dbi(D, F, 0.6), 1e-12);
        // Tiny dish: beamwidth clamped to π.
        let tiny = GaussianBeam::from_dish(1e-3, 1e9, 0.6, 70.0).unwrap();
        assert_eq!(tiny.beamwidth_rad(0.0), PI);
    }

    #[test]
    fn gaussian_rejects_invalid_parameters() {
        assert!(matches!(
            GaussianBeam::new(f64::NAN, 0.1),
            Err(AntennaError::InvalidParameter {
                name: "peak_gain_dbi",
                ..
            })
        ));
        assert!(matches!(
            GaussianBeam::new(30.0, 0.0),
            Err(AntennaError::InvalidParameter {
                name: "bw_major_rad",
                ..
            })
        ));
        assert!(matches!(
            GaussianBeam::elliptical(30.0, 0.1, 4.0, 0.0),
            Err(AntennaError::InvalidParameter {
                name: "bw_minor_rad",
                ..
            })
        ));
        assert!(matches!(
            GaussianBeam::elliptical(30.0, 0.1, 0.1, f64::NAN),
            Err(AntennaError::InvalidParameter {
                name: "orientation_rad",
                ..
            })
        ));
        assert!(GaussianBeam::from_dish(0.0, F, 0.6, 70.0).is_err());
        assert!(GaussianBeam::from_dish(D, 0.0, 0.6, 70.0).is_err());
        assert!(GaussianBeam::from_dish(D, F, 2.0, 70.0).is_err());
        assert!(GaussianBeam::from_dish(D, F, 0.6, 0.0).is_err());
        assert!(GaussianBeam::new(30.0, 0.1).unwrap().with_floor_dbi(f64::NAN).is_err());
    }
}
