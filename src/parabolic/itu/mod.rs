//! ITU-R reference and design-objective antenna patterns.
//!
//! Each formula is transcribed from the text of the cited recommendation
//! (latest edition in force, checked on itu.int). Angles inside the formulas
//! are in degrees, as in the recommendations; the public API takes radians.
//!
//! | Recommendation | Types | Antennas |
//! |---|---|---|
//! | S.465-6 (01/2010) | [`ItuS465`] | FSS earth stations, coordination |
//! | S.580-6 (01/2004) | [`ItuS580`] | FSS earth stations, design objective |
//! | S.1428-1 (02/2001) | [`ItuS1428`] | FSS earth stations, non-GSO interference |
//! | F.699-9 (02/2026) | [`ItuF699`] | fixed wireless systems, side-lobe peaks |
//! | F.1245-3 (01/2019) | [`ItuF1245`] | fixed wireless systems, average side lobes |
//! | S.672-4 (09/1997) | [`ItuS672`] | GSO satellites, single feed |
//! | S.1528-1 (08/2026) | [`ItuS1528`], [`ItuS1528Ngso`], [`ItuS1528Taylor`] | non-GSO satellites |
//! | BO.1213-1 (11/2005) | [`ItuBo1213`] | BSS receiving earth stations |
//! | BO.1295-0 (10/1997) | [`ItuBo1295`] | BSS feeder-link transmitting earth stations |
//! | RR Appendices 30 and 30A (Edition 2024) | [`Ap30EarthReceive`], [`Ap30Satellite`], [`Ap30aEarthTransmit`] | BSS Plans (APS30 / APS30A) |
//!
//! Patterns that the recommendation splits into co-polar and cross-polar
//! components take a [`PolarComponent`]; the others give the co-polar (or
//! single) pattern of the recommendation.

use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite, require_positive};
use crate::math::ellipse_radius;
use crate::parabolic::require_efficiency;

mod bss;
mod earth;
mod fixed;
mod rr_ap30;
mod satellite;

pub use bss::{BO1295_PEAK_GAIN_DBI, ItuBo1213, ItuBo1295};
pub use earth::{ItuS465, ItuS580, ItuS1428, S465Variant};
pub use fixed::{F1245Model, ItuF699, ItuF1245};
pub use rr_ap30::{Ap30EarthKind, Ap30EarthReceive, Ap30Satellite, Ap30SatelliteKind, Ap30aEarthTransmit, Ap30aPlan};
pub use satellite::{ItuS672, ItuS1528, ItuS1528Ngso, ItuS1528Taylor, S672SideLobe, S1528Orbit, S1528SideLobe};

/// Polarisation component of a pattern that the recommendation defines for
/// both the co-polar and the cross-polar response.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PolarComponent {
    /// Co-polar (wanted polarisation) gain.
    #[default]
    CoPolar,
    /// Cross-polar (orthogonal polarisation) gain.
    CrossPolar,
}

/// Main-lobe coefficient of the Radio Regulations Appendix 7/8 pattern,
/// shared by S.465 (crate choice), S.1428, F.699, F.1245 and BO.1213:
/// `G = Gmax − 2.5·10⁻³ (D/λ · φ)²` (φ in degrees).
const MAIN_LOBE_COEFF: f64 = 2.5e-3;

/// Angle (degrees) beyond which the far side-lobe level of S.465, S.580,
/// F.699 and F.1245 (below 70 GHz) applies.
const FAR_ANGLE_DEG: f64 = 48.0;
/// Same angle for F.699 and F.1245 from 70 GHz.
const FAR_ANGLE_HIGH_BAND_DEG: f64 = 120.0;

/// `Gmax − 2.5·10⁻³ (D/λ · φ)²`.
#[inline]
fn main_lobe_parabola_dbi(peak_gain_dbi: f64, d_over_lambda: f64, phi_deg: f64) -> f64 {
    let x = d_over_lambda * phi_deg;
    peak_gain_dbi - MAIN_LOBE_COEFF * x * x
}

/// `φm = (20 λ/D) √(Gmax − G1)` degrees: where the main-lobe parabola meets
/// the first side-lobe level `G1` (S.1428, F.699, F.1245, BO.1213).
fn phi_m_deg(d_over_lambda: f64, peak_gain_dbi: f64, g1_dbi: f64) -> f64 {
    20.0 / d_over_lambda * (peak_gain_dbi - g1_dbi).sqrt()
}

/// Checks `Gmax > G1`, without which `φm` is undefined.
fn require_main_lobe_above(peak_gain_dbi: f64, g1_dbi: f64) -> Result<()> {
    if peak_gain_dbi > g1_dbi {
        Ok(())
    } else {
        Err(invalid(
            "peak_gain_dbi",
            format!("must exceed the first side-lobe gain G1 = {g1_dbi} dBi, got {peak_gain_dbi}"),
        ))
    }
}

/// Main lobe used below `φmin`, where S.465 and S.580 define no gain: the
/// RR Appendix 7 parabola `Gmax − 2.5·10⁻³ (D/λ · φ)²`, not lower than the
/// envelope value at `φmin` and never above `Gmax`. This is a crate choice
/// that keeps the pattern continuous at `φmin`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct MainLobe {
    peak_gain_dbi: f64,
    d_over_lambda: f64,
    phi_min_deg: f64,
    gain_at_phi_min_dbi: f64,
}

impl MainLobe {
    fn gain_dbi(&self, phi_deg: f64) -> f64 {
        main_lobe_parabola_dbi(self.peak_gain_dbi, self.d_over_lambda, phi_deg)
            .max(self.gain_at_phi_min_dbi)
            .min(self.peak_gain_dbi)
    }
}

/// Step envelope shared by F.699 and F.1245 (average model):
///
/// ```text
/// G = Gmax − 2.5·10⁻³ (D/λ · φ)²   φ < φm
///   = G1                           φm ≤ φ < g1_end
///   = near_const − 25 log φ        g1_end ≤ φ < near_end
///   = far                          near_end ≤ φ ≤ 180°
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
struct StepEnvelope {
    peak_gain_dbi: f64,
    d_over_lambda: f64,
    phi_m_deg: f64,
    g1_dbi: f64,
    g1_end_deg: f64,
    near_const_db: f64,
    near_end_deg: f64,
    far_dbi: f64,
}

impl StepEnvelope {
    fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        if phi_deg < self.phi_m_deg {
            main_lobe_parabola_dbi(self.peak_gain_dbi, self.d_over_lambda, phi_deg)
        } else if phi_deg < self.g1_end_deg {
            self.g1_dbi
        } else if phi_deg < self.near_end_deg {
            self.near_const_db - 25.0 * phi_deg.log10()
        } else if phi_deg.is_nan() {
            f64::NAN
        } else {
            self.far_dbi
        }
    }
}

/// Off-boresight angle in degrees, `NaN` for non-finite directions.
#[inline]
fn off_axis_deg(dir: Direction) -> f64 {
    dir.theta_rad().to_degrees()
}

/// Azimuth of the plane containing boresight and `dir`.
#[inline]
fn plane_rad(dir: Direction) -> f64 {
    dir.v().atan2(dir.u())
}

fn require_d_over_lambda(d_over_lambda: f64) -> Result<f64> {
    require_positive("d_over_lambda", d_over_lambda)
}

fn dish_d_over_lambda(diameter_m: f64, freq_hz: f64) -> Result<f64> {
    let d = require_positive("diameter_m", diameter_m)?;
    let f = require_positive("freq_hz", freq_hz)?;
    Ok(d / rf_r::fspl::wavelength_m(f))
}

fn dish_d_over_lambda_and_gain(diameter_m: f64, freq_hz: f64, efficiency: f64) -> Result<(f64, f64)> {
    let dl = dish_d_over_lambda(diameter_m, freq_hz)?;
    let eta = require_efficiency("efficiency", efficiency)?;
    Ok((dl, rf_r::antenna::parabolic_gain_dbi(diameter_m, freq_hz, eta)))
}

/// Elliptical beam cross-section: an angular width (degrees) along the
/// major axis, which lies at azimuth `orientation_rad`, and across it.
/// The width in another plane follows the ellipse
/// ([`ellipse_radius`]). Used for half beamwidths (S.672, S.1528) and full
/// beamwidths (RR Appendix 30).
#[derive(Debug, Clone, Copy, PartialEq)]
struct BeamEllipse {
    major_deg: f64,
    minor_deg: f64,
    orientation_rad: f64,
}

impl BeamEllipse {
    /// Validates `0 < minor ≤ major < max_major_deg` (given in radians) and
    /// a finite orientation. `names` are the parameter names for errors.
    fn new(
        names: (&'static str, &'static str),
        major_rad: f64,
        minor_rad: f64,
        orientation_rad: f64,
        max_major_deg: f64,
    ) -> Result<Self> {
        let major_deg = major_rad.to_degrees();
        let minor_deg = minor_rad.to_degrees();
        if !(major_deg > 0.0 && major_deg < max_major_deg) {
            return Err(invalid(
                names.0,
                format!("must be in (0, {max_major_deg}°), got {major_rad} rad"),
            ));
        }
        if !(minor_deg > 0.0 && minor_deg <= major_deg) {
            return Err(invalid(
                names.1,
                format!("must be in (0, {}], got {minor_rad}", names.0),
            ));
        }
        Ok(Self {
            major_deg,
            minor_deg,
            orientation_rad: require_finite("orientation_rad", orientation_rad)?,
        })
    }

    /// Axis ratio `z = major / minor` (≥ 1).
    fn axis_ratio(&self) -> f64 {
        self.major_deg / self.minor_deg
    }

    /// Width (degrees) in the plane at azimuth `phi_rad`.
    #[inline]
    fn width_deg(&self, phi_rad: f64) -> f64 {
        ellipse_radius(self.major_deg, self.minor_deg, phi_rad - self.orientation_rad)
    }
}

/// `b` of S.672-4 and S.1528-1 Table 1.
const ENVELOPE_B: f64 = 6.32;
/// Far side-lobe level `LF` of S.672-4 and S.1528-1 `recommends 1.2` (dBi).
const ENVELOPE_LF_DBI: f64 = 0.0;

/// Six-segment envelope shared by S.672-4 `recommends 1` and S.1528-1
/// `recommends 1.2` (ψ in degrees):
///
/// ```text
/// G = Gm − 3 (ψ/ψb)^α             ψ ≤ a ψb             (1)
///   = Gm + LN + 20 log z          a ψb < ψ ≤ 0.5 b ψb  (2a)
///   = Gm + LN                     0.5 b ψb < ψ ≤ b ψb  (2b)
///   = X − 25 log ψ                b ψb < ψ ≤ Y         (3)
///   = LF = 0 dBi                  Y < ψ ≤ 90°          (4a)
///   = LB                          90° < ψ ≤ 180°       (4b)
/// X = Gm + LN + 25 log(b ψb),  Y = b ψb 10^(0.04 (Gm + LN − LF))
/// LB = max(15 + LN + 0.25 Gm + 5 log z, 0) dBi
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
struct SideLobeEnvelope {
    peak_gain_dbi: f64,
    half_beamwidth: BeamEllipse,
    level_db: f64,
    a: f64,
    alpha: f64,
    twenty_log_z: f64,
    back_lobe_dbi: f64,
}

impl SideLobeEnvelope {
    /// `a_of_log_z` gives Table 1's `a` from `log z`; it must be positive.
    fn new(
        peak_gain_dbi: f64,
        half_beamwidth: BeamEllipse,
        level_db: f64,
        alpha: f64,
        a_of_log_z: impl Fn(f64) -> f64,
    ) -> Result<Self> {
        let gm = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let log_z = half_beamwidth.axis_ratio().log10();
        let a = a_of_log_z(log_z);
        if a.is_nan() || a <= 0.0 {
            return Err(invalid("half_minor_rad", format!("axis ratio too large: a = {a}")));
        }
        Ok(Self {
            peak_gain_dbi: gm,
            half_beamwidth,
            level_db,
            a,
            alpha,
            twenty_log_z: 20.0 * log_z,
            back_lobe_dbi: (15.0 + level_db + 0.25 * gm + 5.0 * log_z).max(0.0),
        })
    }

    fn gain_at_deg(&self, psi_deg: f64, psi_b_deg: f64) -> f64 {
        let r = psi_deg / psi_b_deg;
        let gm_ln = self.peak_gain_dbi + self.level_db;
        let b_psi = ENVELOPE_B * psi_b_deg;
        if r <= self.a {
            self.peak_gain_dbi - 3.0 * r.powf(self.alpha)
        } else if r <= 0.5 * ENVELOPE_B {
            gm_ln + self.twenty_log_z
        } else if r <= ENVELOPE_B {
            gm_ln
        } else if psi_deg <= b_psi * 10f64.powf(0.04 * (gm_ln - ENVELOPE_LF_DBI)) {
            gm_ln + 25.0 * b_psi.log10() - 25.0 * psi_deg.log10()
        } else if psi_deg <= 90.0 {
            ENVELOPE_LF_DBI
        } else if psi_deg.is_nan() {
            f64::NAN
        } else {
            self.back_lobe_dbi
        }
    }

    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let psi_b = self.half_beamwidth.width_deg(plane_rad(dir));
        self.gain_at_deg(off_axis_deg(dir), psi_b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    #[test]
    fn step_envelope_segments_and_nan() {
        let e = StepEnvelope {
            peak_gain_dbi: 40.0,
            d_over_lambda: 100.0,
            phi_m_deg: 1.0,
            g1_dbi: 20.0,
            g1_end_deg: 2.0,
            near_const_db: 32.0,
            near_end_deg: 48.0,
            far_dbi: -10.0,
        };
        assert_close(e.gain_at_deg(0.5), 40.0 - 2.5e-3 * 2500.0, 1e-12);
        assert_eq!(e.gain_at_deg(1.5), 20.0);
        assert_close(e.gain_at_deg(10.0), 7.0, 1e-12);
        assert_eq!(e.gain_at_deg(48.0), -10.0);
        assert!(e.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn main_lobe_check_and_phi_m() {
        assert!(require_main_lobe_above(30.0, 20.0).is_ok());
        assert!(matches!(
            require_main_lobe_above(20.0, 20.0),
            Err(AntennaError::InvalidParameter {
                name: "peak_gain_dbi",
                ..
            })
        ));
        assert_close(phi_m_deg(100.0, 45.0, 20.0), 1.0, 1e-12);
    }

    #[test]
    fn beam_ellipse_validation() {
        let deg = 1f64.to_radians();
        let e = BeamEllipse::new(("major", "minor"), 2.0 * deg, deg, 0.0, 90.0).unwrap();
        assert_close(e.axis_ratio(), 2.0, 1e-12);
        assert_close(e.width_deg(std::f64::consts::FRAC_PI_2), 1.0, 1e-12);
        assert!(matches!(
            BeamEllipse::new(("major", "minor"), 0.0, deg, 0.0, 90.0),
            Err(AntennaError::InvalidParameter { name: "major", .. })
        ));
        assert!(matches!(
            BeamEllipse::new(("major", "minor"), deg, 2.0 * deg, 0.0, 90.0),
            Err(AntennaError::InvalidParameter { name: "minor", .. })
        ));
        assert!(matches!(
            BeamEllipse::new(("major", "minor"), deg, deg, f64::NAN, 90.0),
            Err(AntennaError::InvalidParameter {
                name: "orientation_rad",
                ..
            })
        ));
    }

    #[test]
    fn dish_helpers_validate() {
        assert!(dish_d_over_lambda(0.0, 1e9).is_err());
        assert!(dish_d_over_lambda(1.0, -1.0).is_err());
        assert!(dish_d_over_lambda_and_gain(1.0, 1e9, 1.5).is_err());
        assert_eq!(PolarComponent::default(), PolarComponent::CoPolar);
    }
}
