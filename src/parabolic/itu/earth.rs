//! FSS earth-station patterns: ITU-R S.465-6, S.580-6 and S.1428-1.

use super::{
    FAR_ANGLE_DEG, MainLobe, dish_d_over_lambda, dish_d_over_lambda_and_gain, main_lobe_parabola_dbi, off_axis_deg,
    phi_m_deg, require_d_over_lambda,
};
use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite};
use crate::pattern::GainPattern;

/// Which pattern of ITU-R S.465-6 applies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum S465Variant {
    /// `recommends 2`: `32 − 25 log φ` from `φmin`, −10 dBi beyond 48°.
    #[default]
    Standard,
    /// Note 5: receiving antennas with `D/λ < 33.3`, `φmin = 2.5°`.
    SmallReceive,
    /// Note 4: antennas with `D/λ ≤ 100` in networks coordinated before
    /// 1993: `52 − 10 log(D/λ) − 25 log φ` from `100 λ/D`,
    /// `10 − 10 log(D/λ)` beyond 48°.
    Pre1993,
}

/// ITU-R S.465-6 reference pattern of FSS earth-station antennas
/// (2–31 GHz), for coordination and interference assessment.
///
/// ```text
/// G = 32 − 25 log φ   dBi   φmin ≤ φ < 48°
///   = −10             dBi   48° ≤ φ ≤ 180°
/// φmin = max(1°, 100 λ/D)               for D/λ ≥ 50
/// φmin = max(2°, 114 (D/λ)^−1.09)       for D/λ < 50
/// ```
///
/// Below `φmin` the crate uses the main lobe described on `MainLobe`
/// (RR Appendix 7 parabola), since S.465 defines none.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuS465 {
    main: MainLobe,
    near_constant_db: f64,
    far_gain_dbi: f64,
}

impl ItuS465 {
    /// Pattern for an antenna of normalised diameter `d_over_lambda` and
    /// peak gain `peak_gain_dbi`.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::{ItuS465, S465Variant}};
    /// let p = ItuS465::new(100.0, 47.7, S465Variant::Standard).unwrap();
    /// let g = p.gain_dbi(Direction::from_theta_phi(10f64.to_radians(), 0.0));
    /// assert!((g - 7.0).abs() < 1e-12);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-positive `D/λ`, a non-finite gain, or a variant whose
    /// condition on `D/λ` is not met (Note 4: `D/λ ≤ 100`; Note 5: `D/λ < 33.3`).
    pub fn new(d_over_lambda: f64, peak_gain_dbi: f64, variant: S465Variant) -> Result<Self> {
        let dl = require_d_over_lambda(d_over_lambda)?;
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let ten_log_dl = 10.0 * dl.log10();
        let (phi_min_deg, near_constant_db, far_gain_dbi) = match variant {
            S465Variant::Standard => (standard_phi_min_deg(dl), 32.0, -10.0),
            S465Variant::SmallReceive => {
                if dl >= 33.3 {
                    return Err(invalid("d_over_lambda", format!("Note 5 needs D/λ < 33.3, got {dl}")));
                }
                (2.5, 32.0, -10.0)
            }
            S465Variant::Pre1993 => {
                if dl > 100.0 {
                    return Err(invalid("d_over_lambda", format!("Note 4 needs D/λ ≤ 100, got {dl}")));
                }
                (100.0 / dl, 52.0 - ten_log_dl, 10.0 - ten_log_dl)
            }
        };
        let side = |phi: f64| near_constant_db - 25.0 * phi.log10();
        Ok(Self {
            main: MainLobe {
                peak_gain_dbi: gmax,
                d_over_lambda: dl,
                phi_min_deg,
                gain_at_phi_min_dbi: side(phi_min_deg),
            },
            near_constant_db,
            far_gain_dbi,
        })
    }

    /// Pattern of a dish: `D/λ` and peak gain from diameter, frequency and
    /// aperture efficiency.
    ///
    /// # Errors
    /// As [`ItuS465::new`], plus invalid dish parameters.
    pub fn from_dish(diameter_m: f64, freq_hz: f64, efficiency: f64, variant: S465Variant) -> Result<Self> {
        let (dl, g) = dish_d_over_lambda_and_gain(diameter_m, freq_hz, efficiency)?;
        Self::new(dl, g, variant)
    }

    /// `φmin` in degrees.
    #[must_use]
    pub fn phi_min_deg(&self) -> f64 {
        self.main.phi_min_deg
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        if phi_deg < self.main.phi_min_deg {
            self.main.gain_dbi(phi_deg)
        } else if phi_deg < FAR_ANGLE_DEG {
            self.near_constant_db - 25.0 * phi_deg.log10()
        } else if phi_deg.is_nan() {
            f64::NAN
        } else {
            self.far_gain_dbi
        }
    }
}

/// `φmin` of S.465-6 `recommends 2` in degrees.
fn standard_phi_min_deg(d_over_lambda: f64) -> f64 {
    if d_over_lambda >= 50.0 {
        (100.0 / d_over_lambda).max(1.0)
    } else {
        (114.0 * d_over_lambda.powf(-1.09)).max(2.0)
    }
}

impl GainPattern for ItuS465 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.main.peak_gain_dbi
    }
}

/// ITU-R S.580-6 design objective for GSO earth-station antennas
/// (`D/λ ≥ 50`), completed with S.465-6 as the recommendation directs:
///
/// ```text
/// G = 29 − 25 log φ   dBi   max(1°, 100 λ/D) ≤ φ ≤ 20°
///   = −3.5            dBi   20° < φ ≤ 26.3°     (Note 5)
///   = 32 − 25 log φ   dBi   26.3° < φ < 48°     (S.465)
///   = −10             dBi   48° ≤ φ ≤ 180°      (S.465)
/// ```
///
/// Below `φmin` the main lobe of [`ItuS465`] is used.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuS580 {
    main: MainLobe,
}

/// End of the 29 − 25 log φ objective (degrees).
const S580_OBJECTIVE_END_DEG: f64 = 20.0;
/// End of the −3.5 dBi bridge to S.465 (degrees).
const S580_BRIDGE_END_DEG: f64 = 26.3;
/// Gain on the bridge between the objective and S.465 (dBi).
const S580_BRIDGE_DBI: f64 = -3.5;

impl ItuS580 {
    /// Design-objective pattern for `D/λ ≥ 50` and peak gain `peak_gain_dbi`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for `D/λ < 50` (Note 3) or a non-finite gain.
    pub fn new(d_over_lambda: f64, peak_gain_dbi: f64) -> Result<Self> {
        let dl = require_d_over_lambda(d_over_lambda)?;
        if dl < 50.0 {
            return Err(invalid("d_over_lambda", format!("S.580 applies to D/λ ≥ 50, got {dl}")));
        }
        let gmax = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let phi_min_deg = (100.0 / dl).max(1.0);
        Ok(Self {
            main: MainLobe {
                peak_gain_dbi: gmax,
                d_over_lambda: dl,
                phi_min_deg,
                gain_at_phi_min_dbi: 29.0 - 25.0 * phi_min_deg.log10(),
            },
        })
    }

    /// Pattern of a dish (see [`ItuS465::from_dish`]).
    ///
    /// # Errors
    /// As [`ItuS580::new`], plus invalid dish parameters.
    pub fn from_dish(diameter_m: f64, freq_hz: f64, efficiency: f64) -> Result<Self> {
        let (dl, g) = dish_d_over_lambda_and_gain(diameter_m, freq_hz, efficiency)?;
        Self::new(dl, g)
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        if phi_deg < self.main.phi_min_deg {
            self.main.gain_dbi(phi_deg)
        } else if phi_deg <= S580_OBJECTIVE_END_DEG {
            29.0 - 25.0 * phi_deg.log10()
        } else if phi_deg <= S580_BRIDGE_END_DEG {
            S580_BRIDGE_DBI
        } else if phi_deg < FAR_ANGLE_DEG {
            32.0 - 25.0 * phi_deg.log10()
        } else if phi_deg.is_nan() {
            f64::NAN
        } else {
            -10.0
        }
    }
}

impl GainPattern for ItuS580 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.main.peak_gain_dbi
    }
}

/// Range of `D/λ` that selects the S.1428-1 pattern.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum S1428Range {
    /// `20 ≤ D/λ ≤ 25`.
    Small,
    /// `25 < D/λ ≤ 100`.
    Medium,
    /// `D/λ > 100`.
    Large,
}

/// Smallest `D/λ` covered by S.1428-1.
const S1428_MIN_D_OVER_LAMBDA: f64 = 20.0;

/// ITU-R S.1428-1 reference FSS earth-station pattern for interference
/// assessment involving non-GSO satellites (10.7–30 GHz). It covers all
/// off-axis angles; φ in degrees, `D` the equivalent diameter.
///
/// ```text
/// 20 ≤ D/λ ≤ 25:   Gmax − 2.5·10⁻³ (D/λ φ)²   0 < φ < φm
///                  G1                         φm ≤ φ < 95 λ/D
///                  29 − 25 log φ              95 λ/D ≤ φ < 33.1°
///                  −9                         33.1° ≤ φ ≤ 80°
///                  −5                         80° < φ ≤ 180°
/// 25 < D/λ ≤ 100:  as above up to 33.1° (included), then
///                  −9 (to 80°), −4 (to 120°), −9 (to 180°)
///   Gmax = 20 log(D/λ) + 7.7,  G1 = 29 − 25 log(95 λ/D)
/// D/λ > 100:       Gmax − 2.5·10⁻³ (D/λ φ)²   0 < φ < φm
///                  G1                         φm ≤ φ < φr
///                  29 − 25 log φ              φr ≤ φ < 10°
///                  34 − 30 log φ              10° ≤ φ < 34.1°
///                  −12, −7, −12               34.1°–80°, 80°–120°, 120°–180°
///   Gmax = 20 log(D/λ) + 8.4,  G1 = −1 + 15 log(D/λ),  φr = 15.85 (D/λ)^−0.6
/// φm = (20 λ/D) √(Gmax − G1)
/// ```
///
/// The text leaves φ = 33.1° undefined for `20 ≤ D/λ ≤ 25`; the crate uses
/// −9 dBi there (the side-lobe formula gives −9.0 dBi too).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuS1428 {
    peak_gain_dbi: f64,
    d_over_lambda: f64,
    phi_m_deg: f64,
    g1_dbi: f64,
    phi_r_deg: f64,
    range: S1428Range,
}

impl ItuS1428 {
    /// Pattern for normalised (equivalent) diameter `d_over_lambda`; the
    /// peak gain follows from `D/λ` as in the recommendation.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern, parabolic::itu::ItuS1428};
    /// let p = ItuS1428::new(50.0).unwrap();
    /// assert!((p.peak_gain_dbi() - (20.0 * 50f64.log10() + 7.7)).abs() < 1e-12);
    /// let g = p.gain_dbi(Direction::from_theta_phi(100f64.to_radians(), 0.0));
    /// assert_eq!(g, -4.0);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for `D/λ < 20` (not covered by the recommendation) or non-finite.
    pub fn new(d_over_lambda: f64) -> Result<Self> {
        let dl = require_d_over_lambda(d_over_lambda)?;
        if !(dl >= S1428_MIN_D_OVER_LAMBDA && dl.is_finite()) {
            return Err(invalid(
                "d_over_lambda",
                format!("S.1428 covers finite D/λ ≥ 20, got {dl}"),
            ));
        }
        let log_dl = dl.log10();
        let (range, peak_gain_dbi, g1_dbi, phi_r_deg) = if dl > 100.0 {
            (
                S1428Range::Large,
                20.0 * log_dl + 8.4,
                -1.0 + 15.0 * log_dl,
                15.85 * dl.powf(-0.6),
            )
        } else {
            let phi_r = 95.0 / dl;
            let range = if dl > 25.0 {
                S1428Range::Medium
            } else {
                S1428Range::Small
            };
            (range, 20.0 * log_dl + 7.7, 29.0 - 25.0 * phi_r.log10(), phi_r)
        };
        Ok(Self {
            peak_gain_dbi,
            d_over_lambda: dl,
            phi_m_deg: phi_m_deg(dl, peak_gain_dbi, g1_dbi),
            g1_dbi,
            phi_r_deg,
            range,
        })
    }

    /// Pattern of a dish of diameter `diameter_m` at `freq_hz`.
    ///
    /// # Errors
    /// As [`ItuS1428::new`], plus a non-positive diameter or frequency.
    pub fn from_dish(diameter_m: f64, freq_hz: f64) -> Result<Self> {
        Self::new(dish_d_over_lambda(diameter_m, freq_hz)?)
    }

    /// `φm` (degrees), end of the main lobe.
    #[must_use]
    pub fn phi_m_deg(&self) -> f64 {
        self.phi_m_deg
    }

    /// First side-lobe gain `G1` (dBi).
    #[must_use]
    pub fn first_side_lobe_dbi(&self) -> f64 {
        self.g1_dbi
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        if phi_deg.is_nan() {
            return f64::NAN;
        }
        if phi_deg < self.phi_m_deg {
            return main_lobe_parabola_dbi(self.peak_gain_dbi, self.d_over_lambda, phi_deg);
        }
        if phi_deg < self.phi_r_deg {
            return self.g1_dbi;
        }
        match self.range {
            S1428Range::Small => small_far_dbi(phi_deg),
            S1428Range::Medium => medium_far_dbi(phi_deg),
            S1428Range::Large => large_far_dbi(phi_deg),
        }
    }
}

/// S.1428-1, `20 ≤ D/λ ≤ 25`, from `95 λ/D`.
fn small_far_dbi(phi_deg: f64) -> f64 {
    if phi_deg < 33.1 {
        29.0 - 25.0 * phi_deg.log10()
    } else if phi_deg <= 80.0 {
        -9.0
    } else {
        -5.0
    }
}

/// S.1428-1, `25 < D/λ ≤ 100`, from `95 λ/D`.
fn medium_far_dbi(phi_deg: f64) -> f64 {
    if phi_deg <= 33.1 {
        29.0 - 25.0 * phi_deg.log10()
    } else if phi_deg <= 80.0 {
        -9.0
    } else if phi_deg <= 120.0 {
        -4.0
    } else {
        -9.0
    }
}

/// S.1428-1, `D/λ > 100`, from `φr`.
fn large_far_dbi(phi_deg: f64) -> f64 {
    if phi_deg < 10.0 {
        29.0 - 25.0 * phi_deg.log10()
    } else if phi_deg < 34.1 {
        34.0 - 30.0 * phi_deg.log10()
    } else if phi_deg < 80.0 {
        -12.0
    } else if phi_deg < 120.0 {
        -7.0
    } else {
        -12.0
    }
}

impl GainPattern for ItuS1428 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
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

    fn at_deg(p: &impl GainPattern, deg: f64) -> f64 {
        p.gain_dbi(Direction::from_theta_phi(deg.to_radians(), 0.3))
    }
    #[test]
    fn s465_phi_min_for_large_and_small_dishes() {
        assert_close(standard_phi_min_deg(200.0), 1.0, 1e-15);
        assert_close(standard_phi_min_deg(80.0), 1.25, 1e-15);
        assert_close(standard_phi_min_deg(50.0), 2.0, 1e-15);
        assert_close(standard_phi_min_deg(49.9), 2.0, 1e-15);
        assert_close(standard_phi_min_deg(20.0), 114.0 * 20f64.powf(-1.09), 1e-12);
    }

    #[test]
    fn s465_side_lobes_and_far_level() {
        let p = ItuS465::new(100.0, 48.0, S465Variant::Standard).unwrap();
        assert_close(at_deg(&p, 1.0), 32.0, 1e-9);
        assert_close(at_deg(&p, 10.0), 7.0, 1e-9);
        assert_close(at_deg(&p, 47.999), 32.0 - 25.0 * 47.999f64.log10(), 1e-9);
        assert_eq!(at_deg(&p, 48.0), -10.0);
        assert_eq!(at_deg(&p, 180.0), -10.0);
        assert_eq!(p.peak_gain_dbi(), 48.0);
        assert_eq!(p.phi_min_deg(), 1.0);
    }

    #[test]
    fn s465_main_lobe_is_continuous_and_bounded() {
        let p = ItuS465::new(100.0, 48.0, S465Variant::Standard).unwrap();
        assert_eq!(at_deg(&p, 0.0), 48.0);
        // Parabola: 48 − 2.5e-3 (100 · 0.5)² = 41.75
        assert_close(at_deg(&p, 0.5), 41.75, 1e-9);
        // Floor at the envelope value at φmin (32 dBi).
        assert_close(at_deg(&p, 0.99), 32.0, 1e-9);
        // Small dish whose envelope exceeds Gmax near φmin: capped at Gmax.
        let small = ItuS465::new(10.0, 10.0, S465Variant::Standard).unwrap();
        assert!(at_deg(&small, 1.0) <= 10.0);
    }

    #[test]
    fn s465_note4_and_note5_variants() {
        let n4 = ItuS465::new(50.0, 40.0, S465Variant::Pre1993).unwrap();
        assert_close(n4.phi_min_deg(), 2.0, 1e-15);
        assert_close(at_deg(&n4, 10.0), 52.0 - 10.0 * 50f64.log10() - 25.0, 1e-9);
        assert_close(at_deg(&n4, 90.0), 10.0 - 10.0 * 50f64.log10(), 1e-9);
        let n5 = ItuS465::new(30.0, 35.0, S465Variant::SmallReceive).unwrap();
        assert_eq!(n5.phi_min_deg(), 2.5);
        assert!(matches!(
            ItuS465::new(33.3, 35.0, S465Variant::SmallReceive),
            Err(AntennaError::InvalidParameter {
                name: "d_over_lambda",
                ..
            })
        ));
        assert!(matches!(
            ItuS465::new(150.0, 35.0, S465Variant::Pre1993),
            Err(AntennaError::InvalidParameter {
                name: "d_over_lambda",
                ..
            })
        ));
        assert_eq!(S465Variant::default(), S465Variant::Standard);
    }

    #[test]
    fn s465_validation_and_nan() {
        assert!(ItuS465::new(0.0, 30.0, S465Variant::Standard).is_err());
        assert!(ItuS465::new(100.0, f64::NAN, S465Variant::Standard).is_err());
        let p = ItuS465::new(100.0, 48.0, S465Variant::Standard).unwrap();
        assert!(p.gain_at_deg(f64::NAN).is_nan());
        assert!(at_deg(&p, f64::NAN).is_nan());
    }

    #[test]
    fn s465_from_dish_derives_d_over_lambda_and_gain() {
        let p = ItuS465::from_dish(2.4, 12.5e9, 0.65, S465Variant::Standard).unwrap();
        let dl = 2.4 / rf_r::fspl::wavelength_m(12.5e9);
        assert_close(p.phi_min_deg(), (100.0 / dl).max(1.0), 1e-12);
        assert_close(
            p.peak_gain_dbi(),
            rf_r::antenna::parabolic_gain_dbi(2.4, 12.5e9, 0.65),
            1e-12,
        );
        assert!(ItuS465::from_dish(0.0, 12.5e9, 0.65, S465Variant::Standard).is_err());
        assert!(ItuS465::from_dish(2.4, 0.0, 0.65, S465Variant::Standard).is_err());
        assert!(ItuS465::from_dish(2.4, 12.5e9, 0.0, S465Variant::Standard).is_err());
    }

    #[test]
    fn s580_segments_and_continuity() {
        let p = ItuS580::new(200.0, 52.0).unwrap();
        assert_close(at_deg(&p, 1.0), 29.0, 1e-9);
        assert_close(at_deg(&p, 20.0), 29.0 - 25.0 * 20f64.log10(), 1e-9);
        assert_eq!(at_deg(&p, 20.5), -3.5);
        assert_eq!(at_deg(&p, 26.3), -3.5);
        assert_close(at_deg(&p, 26.31), 32.0 - 25.0 * 26.31f64.log10(), 1e-9);
        // Bridge level meets S.465 at 26.3° (within rounding of the text).
        assert_close(32.0 - 25.0 * 26.3f64.log10(), -3.5, 0.01);
        assert_eq!(at_deg(&p, 60.0), -10.0);
        assert_eq!(at_deg(&p, 0.0), 52.0);
        assert_eq!(p.peak_gain_dbi(), 52.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn s580_validation_and_from_dish() {
        assert!(matches!(
            ItuS580::new(49.0, 40.0),
            Err(AntennaError::InvalidParameter {
                name: "d_over_lambda",
                ..
            })
        ));
        assert!(ItuS580::new(-1.0, 40.0).is_err());
        assert!(ItuS580::new(60.0, f64::INFINITY).is_err());
        let p = ItuS580::from_dish(2.4, 14e9, 0.7).unwrap();
        assert_close(
            p.peak_gain_dbi(),
            rf_r::antenna::parabolic_gain_dbi(2.4, 14e9, 0.7),
            1e-12,
        );
        assert!(ItuS580::from_dish(0.5, 4e9, 0.7).is_err());
        assert!(ItuS580::from_dish(f64::NAN, 4e9, 0.7).is_err());
    }

    #[test]
    fn s1428_small_dish_segments() {
        let dl = 22.0;
        let p = ItuS1428::new(dl).unwrap();
        let gmax = 20.0 * dl.log10() + 7.7;
        let phi_r = 95.0 / dl;
        let g1 = 29.0 - 25.0 * phi_r.log10();
        assert_close(p.peak_gain_dbi(), gmax, 1e-12);
        assert_close(p.first_side_lobe_dbi(), g1, 1e-12);
        assert_close(p.phi_m_deg(), 20.0 / dl * (gmax - g1).sqrt(), 1e-12);
        assert_eq!(at_deg(&p, 0.0), gmax);
        // Main lobe meets G1 at φm.
        assert_close(p.gain_at_deg(p.phi_m_deg() * (1.0 - 1e-12)), g1, 1e-6);
        assert_close(p.gain_at_deg(phi_r), 29.0 - 25.0 * phi_r.log10(), 1e-12);
        assert_close(p.gain_at_deg(20.0), 29.0 - 25.0 * 20f64.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(33.1), -9.0);
        assert_eq!(p.gain_at_deg(80.0), -9.0);
        assert_eq!(p.gain_at_deg(80.1), -5.0);
        assert_eq!(at_deg(&p, 180.0), -5.0);
    }

    #[test]
    fn s1428_medium_dish_back_lobes() {
        let p = ItuS1428::new(60.0).unwrap();
        assert_close(p.gain_at_deg(33.1), 29.0 - 25.0 * 33.1f64.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(33.2), -9.0);
        assert_eq!(p.gain_at_deg(80.0), -9.0);
        assert_eq!(p.gain_at_deg(100.0), -4.0);
        assert_eq!(p.gain_at_deg(120.0), -4.0);
        assert_eq!(p.gain_at_deg(150.0), -9.0);
        // D/λ = 25 is still in the small range: −5 dBi back lobe.
        assert_eq!(ItuS1428::new(25.0).unwrap().gain_at_deg(100.0), -5.0);
    }

    #[test]
    fn s1428_large_dish_segments() {
        let dl = 200.0;
        let p = ItuS1428::new(dl).unwrap();
        let gmax = 20.0 * dl.log10() + 8.4;
        let g1 = -1.0 + 15.0 * dl.log10();
        let phi_r = 15.85 * dl.powf(-0.6);
        assert_close(p.peak_gain_dbi(), gmax, 1e-12);
        assert_close(p.first_side_lobe_dbi(), g1, 1e-12);
        assert!(p.phi_m_deg() < phi_r);
        assert_eq!(p.gain_at_deg(f64::midpoint(p.phi_m_deg(), phi_r)), g1);
        assert_close(p.gain_at_deg(phi_r), 29.0 - 25.0 * phi_r.log10(), 1e-12);
        assert_close(p.gain_at_deg(10.0), 4.0, 1e-12);
        assert_close(p.gain_at_deg(9.999_999), 4.0, 1e-5);
        // 34 − 30 log 34.1 ≈ −12.0: continuous with the −12 dBi plateau.
        assert_close(p.gain_at_deg(34.0), 34.0 - 30.0 * 34f64.log10(), 1e-12);
        assert_close(34.0 - 30.0 * 34.1f64.log10(), -12.0, 0.02);
        assert_eq!(p.gain_at_deg(34.1), -12.0);
        assert_eq!(p.gain_at_deg(80.0), -7.0);
        assert_eq!(p.gain_at_deg(120.0), -12.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn s1428_validation_and_from_dish() {
        for bad in [19.9, 0.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                ItuS1428::new(bad),
                Err(AntennaError::InvalidParameter {
                    name: "d_over_lambda",
                    ..
                })
            ));
        }
        let p = ItuS1428::from_dish(0.6, 20e9).unwrap();
        let dl = 0.6 / rf_r::fspl::wavelength_m(20e9);
        assert_close(p.peak_gain_dbi(), 20.0 * dl.log10() + 7.7, 1e-12);
        assert!(ItuS1428::from_dish(0.0, 20e9).is_err());
    }
}
