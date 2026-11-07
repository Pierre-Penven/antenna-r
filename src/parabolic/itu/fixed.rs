//! Fixed wireless system (point-to-point) patterns: ITU-R F.699-9 (side-lobe
//! peaks) and F.1245-3 (average side lobes).

use super::{
    FAR_ANGLE_DEG, FAR_ANGLE_HIGH_BAND_DEG, StepEnvelope, main_lobe_parabola_dbi, off_axis_deg, phi_m_deg,
    require_d_over_lambda, require_main_lobe_above,
};
use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite, require_positive};
use crate::pattern::GainPattern;

/// `D/λ` above which the "large antenna" equations apply (F.699, F.1245).
const LARGE_D_OVER_LAMBDA: f64 = 100.0;
/// Start of the 70–86 GHz equations (F.699 also uses them up to 174.8 GHz).
const HIGH_BAND_HZ: f64 = 70e9;
/// Start of the 1–70 GHz equations.
const MICROWAVE_HZ: f64 = 1e9;

/// `Gmax ≈ 20 log(D/λ) + 7.7` dBi (F.699-9 `recommends 3`, F.1245-3 Note 2).
fn gain_from_d_over_lambda(d_over_lambda: f64) -> f64 {
    20.0 * d_over_lambda.log10() + 7.7
}

/// `D/λ ≈ 10^((Gmax − 7.7)/20)` (same relation, inverted).
fn d_over_lambda_from_gain(peak_gain_dbi: f64) -> f64 {
    10f64.powf((peak_gain_dbi - 7.7) / 20.0)
}

/// First side-lobe gain `G1 = 2 + 15 log(D/λ)` dBi (F.699, F.1245).
fn first_side_lobe_dbi(d_over_lambda: f64) -> f64 {
    2.0 + 15.0 * d_over_lambda.log10()
}

fn require_frequency(freq_hz: f64, min_hz: f64, max_hz: f64, rec: &str) -> Result<f64> {
    let f = require_positive("freq_hz", freq_hz)?;
    if (min_hz..=max_hz).contains(&f) {
        Ok(f)
    } else {
        Err(invalid(
            "freq_hz",
            format!("{rec} covers {min_hz} Hz to {max_hz} Hz, got {f}"),
        ))
    }
}

/// Lowest frequency of F.699-9 (100 MHz).
const F699_MIN_HZ: f64 = 100e6;
/// Highest frequency of F.699-9 (174.8 GHz).
const F699_MAX_HZ: f64 = 174.8e9;
/// Smallest `D/λ` of F.699-9 `recommends 2.3` (100 MHz–1 GHz).
const F699_UHF_MIN_D_OVER_LAMBDA: f64 = 0.63;

/// ITU-R F.699-9 reference pattern of fixed wireless system antennas
/// (side-lobe peaks, 100 MHz–174.8 GHz), φ in degrees:
///
/// ```text
/// all cases:  Gmax − 2.5·10⁻³ (D/λ φ)²   0 < φ < φm
/// D/λ > 100, 1–70 GHz (2.1.1):   G1 to φr; 32 − 25 log φ to 48°; −10 to 180°
/// D/λ > 100, 70–86 GHz (2.1.2):  G1 to φr; 32 − 25 log φ to 120°; −20 to 180°
/// D/λ ≤ 100, 1–70 GHz (2.2.1):   G1 to 100 λ/D; 52 − 10 log(D/λ) − 25 log φ
///                                to 48°; 10 − 10 log(D/λ) to 180°
/// D/λ ≤ 100, 70–86 GHz (2.2.2):  G1 to 100 λ/D; 52 − 10 log(D/λ) − 25 log φ
///                                to 120°; −10 log(D/λ) to 180°
/// 100 MHz–1 GHz, D/λ > 0.63 (2.3): G1 to 100 λ/D; 52 − 10 log(D/λ) − 25 log φ
///                                to φs; −2 − 5 log(D/λ) to 180°
/// G1 = 2 + 15 log(D/λ),  φm = (20 λ/D) √(Gmax − G1),
/// φr = 15.85 (D/λ)^−0.6,  φs = 144.5 (D/λ)^−0.2
/// ```
///
/// From 86 to 174.8 GHz the 70–86 GHz equations apply (`recommends 2.4`).
/// The pattern is for one (co-)polarisation (Note 4).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuF699 {
    envelope: StepEnvelope,
}

impl ItuF699 {
    /// Pattern for `D/λ`, peak gain and frequency.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::ItuF699};
    /// let p = ItuF699::from_d_over_lambda(200.0, 18e9).unwrap();
    /// let g = p.gain_dbi(Direction::from_theta_phi(10f64.to_radians(), 0.0));
    /// assert!((g - 7.0).abs() < 1e-12); // 32 − 25 log 10
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a frequency outside 100 MHz–174.8 GHz, a non-positive `D/λ`,
    /// `D/λ ≤ 0.63` below 1 GHz, or a peak gain not above `G1`.
    pub fn new(d_over_lambda: f64, peak_gain_dbi: f64, freq_hz: f64) -> Result<Self> {
        let f = require_frequency(freq_hz, F699_MIN_HZ, F699_MAX_HZ, "F.699")?;
        let dl = require_d_over_lambda(d_over_lambda)?;
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let g1 = first_side_lobe_dbi(dl);
        require_main_lobe_above(gmax, g1)?;
        let ten_log_dl = 10.0 * dl.log10();
        let (g1_end_deg, near_const_db, near_end_deg, far_dbi) = if f < MICROWAVE_HZ {
            if dl <= F699_UHF_MIN_D_OVER_LAMBDA {
                return Err(invalid(
                    "d_over_lambda",
                    format!("F.699 below 1 GHz needs D/λ > 0.63, got {dl}"),
                ));
            }
            (
                100.0 / dl,
                52.0 - ten_log_dl,
                144.5 * dl.powf(-0.2),
                -2.0 - 0.5 * ten_log_dl,
            )
        } else {
            let high = f >= HIGH_BAND_HZ;
            let near_end = if high { FAR_ANGLE_HIGH_BAND_DEG } else { FAR_ANGLE_DEG };
            match (dl > LARGE_D_OVER_LAMBDA, high) {
                (true, false) => (15.85 * dl.powf(-0.6), 32.0, near_end, -10.0),
                (true, true) => (15.85 * dl.powf(-0.6), 32.0, near_end, -20.0),
                (false, false) => (100.0 / dl, 52.0 - ten_log_dl, near_end, 10.0 - ten_log_dl),
                (false, true) => (100.0 / dl, 52.0 - ten_log_dl, near_end, -ten_log_dl),
            }
        };
        Ok(Self {
            envelope: StepEnvelope {
                peak_gain_dbi: gmax,
                d_over_lambda: dl,
                phi_m_deg: phi_m_deg(dl, gmax, g1),
                g1_dbi: g1,
                g1_end_deg,
                near_const_db,
                near_end_deg,
                far_dbi,
            },
        })
    }

    /// Pattern when only `D/λ` is known: `Gmax = 20 log(D/λ) + 7.7`
    /// (`recommends 3`).
    ///
    /// # Errors
    /// As [`ItuF699::new`].
    pub fn from_d_over_lambda(d_over_lambda: f64, freq_hz: f64) -> Result<Self> {
        let dl = require_d_over_lambda(d_over_lambda)?;
        Self::new(dl, gain_from_d_over_lambda(dl), freq_hz)
    }

    /// Pattern when only the peak gain is known: `20 log(D/λ) = Gmax − 7.7`
    /// (`recommends 3`).
    ///
    /// # Errors
    /// As [`ItuF699::new`].
    pub fn from_peak_gain_dbi(peak_gain_dbi: f64, freq_hz: f64) -> Result<Self> {
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        Self::new(d_over_lambda_from_gain(gmax), gmax, freq_hz)
    }

    /// `φm` (degrees), end of the main lobe.
    #[must_use]
    pub fn phi_m_deg(&self) -> f64 {
        self.envelope.phi_m_deg
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        self.envelope.gain_at_deg(phi_deg)
    }
}

impl GainPattern for ItuF699 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.envelope.peak_gain_dbi
    }
}

/// Lowest frequency of F.1245-3 (1 GHz).
const F1245_MIN_HZ: f64 = 1e9;
/// Highest frequency of F.1245-3 (86 GHz).
const F1245_MAX_HZ: f64 = 86e9;

/// Which F.1245-3 model to evaluate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum F1245Model {
    /// `recommends 2`: average side-lobe levels, for aggregate interference.
    #[default]
    Average,
    /// Annex 1 (provisional): generalised pattern whose side lobes vary
    /// sinusoidally around the average, for spatial statistical analyses
    /// with a few interference entries only.
    Generalized,
}

/// ITU-R F.1245-3 mathematical model of average (and generalised) radiation
/// patterns of point-to-point fixed wireless antennas, 1–86 GHz, φ in
/// degrees.
///
/// Average model (`recommends 2`):
///
/// ```text
/// all cases:  Gmax − 2.5·10⁻³ (D/λ φ)²   0 < φ < φm
/// D/λ > 100:  G1 to max(φm, φr); 29 − 25 log φ to 48° (120° from 70 GHz);
///             −13 dBi (−23 dBi from 70 GHz) to 180°
/// D/λ ≤ 100:  39 − 5 log(D/λ) − 25 log φ to 48° (120°);
///             −3 − 5 log(D/λ) (−13 − 5 log(D/λ)) to 180°
/// G1 = 2 + 15 log(D/λ),  φm = (20 λ/D) √(Gmax − G1),  φr = 12.02 (D/λ)^−0.6
/// ```
///
/// Generalised model (Annex 1), with
/// `F(φ) = 10 log(0.9 sin²(3πφ / (2φr)) + 0.1)` (≤ 0, 0 at side-lobe peaks):
///
/// ```text
/// max[Gmax − 2.5·10⁻³ (D/λ φ)², G1 + F(φ)]    0 ≤ φ < φr
/// D/λ > 100:  32 − 25 log φ + F(φ) to 48° (120°); −10 + F(φ) (−20 + F(φ)) to 180°
///             φr = 15.85 (D/λ)^−0.6
/// D/λ ≤ 100:  42 − 5 log(D/λ) − 25 log φ + F(φ) to 48° (120°);
///             −5 log(D/λ) + F(φ) (−10 − 5 log(D/λ) + F(φ)) to 180°
///             φr = 39.8 (D/λ)^−0.8
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuF1245 {
    model: F1245Model,
    envelope: StepEnvelope,
    phi_r_deg: f64,
}

impl ItuF1245 {
    /// Pattern for `D/λ`, peak gain, frequency and model.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::{F1245Model, ItuF1245}};
    /// let p = ItuF1245::from_d_over_lambda(200.0, 18e9, F1245Model::Average).unwrap();
    /// let g = p.gain_dbi(Direction::from_theta_phi(10f64.to_radians(), 0.0));
    /// assert!((g - 4.0).abs() < 1e-12); // 29 − 25 log 10
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a frequency outside 1–86 GHz, a non-positive `D/λ` or a peak gain
    /// not above `G1`.
    pub fn new(d_over_lambda: f64, peak_gain_dbi: f64, freq_hz: f64, model: F1245Model) -> Result<Self> {
        let f = require_frequency(freq_hz, F1245_MIN_HZ, F1245_MAX_HZ, "F.1245")?;
        let dl = require_d_over_lambda(d_over_lambda)?;
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let g1 = first_side_lobe_dbi(dl);
        require_main_lobe_above(gmax, g1)?;
        let phi_m = phi_m_deg(dl, gmax, g1);
        let high = f >= HIGH_BAND_HZ;
        let near_end_deg = if high { FAR_ANGLE_HIGH_BAND_DEG } else { FAR_ANGLE_DEG };
        let five_log_dl = 5.0 * dl.log10();
        let large = dl > LARGE_D_OVER_LAMBDA;
        let (g1_end_deg, near_const_db, far_dbi, phi_r_deg) = match (model, large) {
            (F1245Model::Average, true) => (
                phi_m.max(12.02 * dl.powf(-0.6)),
                29.0,
                if high { -23.0 } else { -13.0 },
                f64::NAN,
            ),
            (F1245Model::Average, false) => (
                phi_m,
                39.0 - five_log_dl,
                if high { -13.0 } else { -3.0 } - five_log_dl,
                f64::NAN,
            ),
            (F1245Model::Generalized, true) => {
                (f64::NAN, 32.0, if high { -20.0 } else { -10.0 }, 15.85 * dl.powf(-0.6))
            }
            (F1245Model::Generalized, false) => (
                f64::NAN,
                42.0 - five_log_dl,
                if high { -10.0 } else { 0.0 } - five_log_dl,
                39.8 * dl.powf(-0.8),
            ),
        };
        Ok(Self {
            model,
            envelope: StepEnvelope {
                peak_gain_dbi: gmax,
                d_over_lambda: dl,
                phi_m_deg: phi_m,
                g1_dbi: g1,
                g1_end_deg,
                near_const_db,
                near_end_deg,
                far_dbi,
            },
            phi_r_deg,
        })
    }

    /// Pattern when only `D/λ` is known: `Gmax = 20 log(D/λ) + 7.7` (Note 2).
    ///
    /// # Errors
    /// As [`ItuF1245::new`].
    pub fn from_d_over_lambda(d_over_lambda: f64, freq_hz: f64, model: F1245Model) -> Result<Self> {
        let dl = require_d_over_lambda(d_over_lambda)?;
        Self::new(dl, gain_from_d_over_lambda(dl), freq_hz, model)
    }

    /// Pattern when only the peak gain is known (Note 2, inverted).
    ///
    /// # Errors
    /// As [`ItuF1245::new`].
    pub fn from_peak_gain_dbi(peak_gain_dbi: f64, freq_hz: f64, model: F1245Model) -> Result<Self> {
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        Self::new(d_over_lambda_from_gain(gmax), gmax, freq_hz, model)
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        match self.model {
            F1245Model::Average => self.envelope.gain_at_deg(phi_deg),
            F1245Model::Generalized => self.generalized_at_deg(phi_deg),
        }
    }

    /// Annex 1 generalised pattern.
    fn generalized_at_deg(&self, phi_deg: f64) -> f64 {
        let e = &self.envelope;
        let s = (1.5 * std::f64::consts::PI * phi_deg / self.phi_r_deg).sin();
        let ripple_db = 10.0 * 0.9f64.mul_add(s * s, 0.1).log10();
        if phi_deg < self.phi_r_deg {
            main_lobe_parabola_dbi(e.peak_gain_dbi, e.d_over_lambda, phi_deg).max(e.g1_dbi + ripple_db)
        } else if phi_deg < e.near_end_deg {
            e.near_const_db - 25.0 * phi_deg.log10() + ripple_db
        } else {
            // NaN propagates through `ripple_db`.
            e.far_dbi + ripple_db
        }
    }
}

impl GainPattern for ItuF1245 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.envelope.peak_gain_dbi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    fn at_deg(p: &impl GainPattern, deg: f64) -> f64 {
        p.gain_dbi(Direction::from_theta_phi(deg.to_radians(), 1.0))
    }

    #[test]
    fn f699_large_antenna_low_band() {
        let dl = 200.0_f64;
        let p = ItuF699::from_d_over_lambda(dl, 18e9).unwrap();
        let gmax = 20.0 * dl.log10() + 7.7;
        let g1 = 2.0 + 15.0 * dl.log10();
        let phi_r = 15.85 * dl.powf(-0.6);
        assert_close(p.peak_gain_dbi(), gmax, 1e-12);
        assert_close(p.phi_m_deg(), 20.0 / dl * (gmax - g1).sqrt(), 1e-12);
        assert_eq!(at_deg(&p, 0.0), gmax);
        assert_close(p.gain_at_deg(0.1), gmax - 2.5e-3 * 400.0, 1e-12);
        assert_eq!(p.gain_at_deg(f64::midpoint(p.phi_m_deg(), phi_r)), g1);
        assert_close(p.gain_at_deg(phi_r), 32.0 - 25.0 * phi_r.log10(), 1e-12);
        assert_close(p.gain_at_deg(47.9), 32.0 - 25.0 * 47.9f64.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(48.0), -10.0);
        assert_eq!(p.gain_at_deg(180.0), -10.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn f699_large_antenna_high_band_and_above_86_ghz() {
        for f in [70e9, 80e9, 150e9, 174.8e9] {
            let p = ItuF699::from_d_over_lambda(300.0, f).unwrap();
            assert_close(p.gain_at_deg(100.0), -18.0, 1e-12);
            assert_eq!(p.gain_at_deg(120.0), -20.0);
        }
    }

    #[test]
    fn f699_small_antenna_both_bands_are_continuous() {
        let dl = 50.0_f64;
        let x = dl.log10();
        let low = ItuF699::from_d_over_lambda(dl, 10e9).unwrap();
        assert_eq!(low.gain_at_deg(1.99), 2.0 + 15.0 * x);
        assert_close(low.gain_at_deg(2.0), 52.0 - 10.0 * x - 25.0 * 2f64.log10(), 1e-12);
        assert_close(low.gain_at_deg(48.0), 10.0 - 10.0 * x, 1e-12);
        // Continuity at 48°: 52 − 25 log 48 ≈ 10.
        assert_close(low.gain_at_deg(47.999_999), low.gain_at_deg(48.0), 0.05);
        let high = ItuF699::from_d_over_lambda(dl, 75e9).unwrap();
        assert_close(high.gain_at_deg(100.0), 52.0 - 10.0 * x - 50.0, 1e-12);
        assert_close(high.gain_at_deg(120.0), -10.0 * x, 1e-12);
        assert_close(high.gain_at_deg(119.999_999), high.gain_at_deg(120.0), 0.05);
    }

    #[test]
    fn f699_uhf_band() {
        let dl = 5.0_f64;
        let x = dl.log10();
        let p = ItuF699::from_d_over_lambda(dl, 400e6).unwrap();
        let phi_s = 144.5 * dl.powf(-0.2);
        assert_close(p.gain_at_deg(30.0), 52.0 - 10.0 * x - 25.0 * 30f64.log10(), 1e-12);
        assert_close(p.gain_at_deg(phi_s), -2.0 - 5.0 * x, 1e-12);
        // −2 − 5 log(D/λ) is the side-lobe formula at φs (to the rounding of 144.5).
        assert_close(52.0 - 10.0 * x - 25.0 * phi_s.log10(), -2.0 - 5.0 * x, 0.01);
        assert!(matches!(
            ItuF699::from_d_over_lambda(0.63, 400e6),
            Err(AntennaError::InvalidParameter {
                name: "d_over_lambda",
                ..
            })
        ));
    }

    #[test]
    fn f699_constructors_and_validation() {
        let p = ItuF699::from_peak_gain_dbi(48.0, 20e9).unwrap();
        assert_eq!(p.peak_gain_dbi(), 48.0);
        let q = ItuF699::from_d_over_lambda(10f64.powf((48.0 - 7.7) / 20.0), 20e9).unwrap();
        assert_close(p.gain_at_deg(3.0), q.gain_at_deg(3.0), 1e-9);
        for f in [99e6, 175e9, 0.0, f64::NAN] {
            assert!(matches!(
                ItuF699::new(100.0, 48.0, f),
                Err(AntennaError::InvalidParameter { name: "freq_hz", .. })
            ));
        }
        assert!(matches!(
            ItuF699::new(100.0, 20.0, 10e9),
            Err(AntennaError::InvalidParameter {
                name: "peak_gain_dbi",
                ..
            })
        ));
        assert!(ItuF699::new(0.0, 48.0, 10e9).is_err());
        assert!(ItuF699::new(100.0, f64::NAN, 10e9).is_err());
        assert!(ItuF699::from_d_over_lambda(-1.0, 10e9).is_err());
        assert!(ItuF699::from_peak_gain_dbi(f64::INFINITY, 10e9).is_err());
    }

    #[test]
    fn f1245_average_large_antenna() {
        let dl = 200.0_f64;
        let p = ItuF1245::from_d_over_lambda(dl, 18e9, F1245Model::Average).unwrap();
        let g1 = 2.0 + 15.0 * dl.log10();
        let phi_r = 12.02 * dl.powf(-0.6);
        // φm ≈ 0.415° < 0.45° < φr ≈ 0.501°: first side-lobe plateau.
        assert_eq!(p.gain_at_deg(0.45), g1);
        assert_close(p.gain_at_deg(phi_r), 29.0 - 25.0 * phi_r.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(48.0), -13.0);
        let high = ItuF1245::from_d_over_lambda(dl, 80e9, F1245Model::Average).unwrap();
        assert_close(high.gain_at_deg(100.0), -21.0, 1e-12);
        assert_eq!(high.gain_at_deg(120.0), -23.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
        assert_eq!(at_deg(&p, 0.0), p.peak_gain_dbi());
    }

    #[test]
    fn f1245_average_is_3_db_below_f699_peaks() {
        // Note 1: § 2.1 average side lobes are 3 dB below F.699 § 2.1.
        let avg = ItuF1245::from_d_over_lambda(300.0, 18e9, F1245Model::Average).unwrap();
        let peak = ItuF699::from_d_over_lambda(300.0, 18e9).unwrap();
        assert_close(peak.gain_at_deg(20.0) - avg.gain_at_deg(20.0), 3.0, 1e-12);
    }

    #[test]
    fn f1245_average_small_antenna() {
        let dl = 40.0_f64;
        let x = dl.log10();
        let p = ItuF1245::from_d_over_lambda(dl, 10e9, F1245Model::Average).unwrap();
        let phi_m = 20.0 / dl * (20.0 * x + 7.7 - 2.0 - 15.0 * x).sqrt();
        assert_close(p.gain_at_deg(phi_m), 39.0 - 5.0 * x - 25.0 * phi_m.log10(), 1e-12);
        assert_close(p.gain_at_deg(48.0), -3.0 - 5.0 * x, 1e-12);
        let high = ItuF1245::from_d_over_lambda(dl, 80e9, F1245Model::Average).unwrap();
        assert_close(high.gain_at_deg(120.0), -13.0 - 5.0 * x, 1e-12);
    }

    #[test]
    fn f1245_generalized_large_antenna() {
        let dl = 200.0_f64;
        let p = ItuF1245::from_d_over_lambda(dl, 18e9, F1245Model::Generalized).unwrap();
        let phi_r = 15.85 * dl.powf(-0.6);
        let g1 = 2.0 + 15.0 * dl.log10();
        let ripple = |phi: f64| 10.0 * (0.9 * (1.5 * std::f64::consts::PI * phi / phi_r).sin().powi(2) + 0.1).log10();
        assert_eq!(p.gain_at_deg(0.0), p.peak_gain_dbi());
        // At φr the phase is 1.5π: a side-lobe peak, F = 0.
        assert_close(ripple(phi_r), 0.0, 1e-12);
        assert_close(p.gain_at_deg(phi_r), 32.0 - 25.0 * phi_r.log10(), 1e-12);
        let near = 0.9 * phi_r;
        let expected = (p.peak_gain_dbi() - 2.5e-3 * (dl * near).powi(2)).max(g1 + ripple(near));
        assert_close(p.gain_at_deg(near), expected, 1e-12);
        assert_close(p.gain_at_deg(10.0), 7.0 + ripple(10.0), 1e-12);
        assert_close(p.gain_at_deg(60.0), -10.0 + ripple(60.0), 1e-12);
        // F never falls below −10 dB.
        assert!(ripple(2.0 * phi_r / 3.0) >= -10.0 - 1e-12);
        let high = ItuF1245::from_d_over_lambda(dl, 80e9, F1245Model::Generalized).unwrap();
        assert_close(high.gain_at_deg(150.0), -20.0 + ripple(150.0), 1e-12);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn f1245_generalized_small_antenna() {
        let dl = 50.0_f64;
        let x = dl.log10();
        let phi_r = 39.8 * dl.powf(-0.8);
        let ripple = |phi: f64| 10.0 * (0.9 * (1.5 * std::f64::consts::PI * phi / phi_r).sin().powi(2) + 0.1).log10();
        let low = ItuF1245::from_d_over_lambda(dl, 10e9, F1245Model::Generalized).unwrap();
        assert_close(
            low.gain_at_deg(20.0),
            42.0 - 5.0 * x - 25.0 * 20f64.log10() + ripple(20.0),
            1e-12,
        );
        assert_close(low.gain_at_deg(90.0), -5.0 * x + ripple(90.0), 1e-12);
        let high = ItuF1245::from_d_over_lambda(dl, 80e9, F1245Model::Generalized).unwrap();
        assert_close(high.gain_at_deg(150.0), -10.0 - 5.0 * x + ripple(150.0), 1e-12);
    }

    #[test]
    fn f1245_constructors_and_validation() {
        assert_eq!(F1245Model::default(), F1245Model::Average);
        let p = ItuF1245::from_peak_gain_dbi(45.0, 20e9, F1245Model::Average).unwrap();
        assert_eq!(p.peak_gain_dbi(), 45.0);
        for f in [0.9e9, 87e9] {
            assert!(matches!(
                ItuF1245::new(100.0, 45.0, f, F1245Model::Average),
                Err(AntennaError::InvalidParameter { name: "freq_hz", .. })
            ));
        }
        assert!(ItuF1245::new(100.0, 10.0, 10e9, F1245Model::Average).is_err());
        assert!(ItuF1245::new(-1.0, 45.0, 10e9, F1245Model::Average).is_err());
        assert!(ItuF1245::new(100.0, f64::NAN, 10e9, F1245Model::Average).is_err());
        assert!(ItuF1245::from_d_over_lambda(0.0, 10e9, F1245Model::Average).is_err());
        assert!(ItuF1245::from_peak_gain_dbi(f64::NAN, 10e9, F1245Model::Average).is_err());
    }
}
