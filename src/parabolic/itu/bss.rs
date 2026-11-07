//! Broadcasting-satellite service earth-station patterns: ITU-R BO.1213-1
//! (receiving, 11.7–12.75 GHz) and BO.1295-0 (feeder-link transmitting,
//! off-axis e.i.r.p.).

use super::{
    PolarComponent, dish_d_over_lambda, main_lobe_parabola_dbi, off_axis_deg, phi_m_deg, require_d_over_lambda,
};
use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite};
use crate::parabolic::require_efficiency;
use crate::pattern::GainPattern;

/// Smallest `D/λ` covered by BO.1213-1.
const BO1213_MIN_D_OVER_LAMBDA: f64 = 11.0;
/// Start of the −5 dBi plateau of the co-polar pattern, `10^(34/25)` degrees.
const BO1213_PHI_B_DEG: f64 = 22.908_676_527_677_73;
/// Start of the −5 dBi plateau of the cross-polar pattern, `10^(26/25)` degrees.
const BO1213_PHI_2_DEG: f64 = 10.964_781_961_431_85;
/// Start of the 0 dBi back region (degrees).
const BO1213_BACK_DEG: f64 = 70.0;

/// ITU-R BO.1213-1 reference receiving earth-station pattern for the BSS
/// (11.7–12.75 GHz), valid for `D/λ ≥ 11`, φ in degrees:
///
/// ```text
/// Co-polar:
///   Gmax − 2.5·10⁻³ (D/λ φ)²     0 ≤ φ < φm
///   G1 = 29 − 25 log φr          φm ≤ φ < φr,   φr = 95 λ/D
///   29 − 25 log φ                φr ≤ φ < φb,   φb = 10^(34/25)
///   −5 dBi                       φb ≤ φ < 70°
///   0 dBi                        70° ≤ φ ≤ 180°
/// Cross-polar:
///   Gmax − 25                                0 ≤ φ < 0.25 φ0
///   Gmax − 25 + 8 (φ − 0.25 φ0)/(0.19 φ0)    0.25 φ0 ≤ φ < 0.44 φ0
///   Gmax − 17                                0.44 φ0 ≤ φ < φ0
///   Gmax − 17 + C (φ − φ0)/(φ1 − φ0)         φ0 ≤ φ < φ1
///   21 − 25 log φ                            φ1 ≤ φ < φ2,   φ2 = 10^(26/25)
///   −5 dBi                                   φ2 ≤ φ < 70°
///   0 dBi                                    70° ≤ φ ≤ 180°
/// Gmax = 10 log(η (π D/λ)²),  φm = (λ/D) √((Gmax − G1)/0.0025)
/// φ0 = 2 (λ/D) √(3/0.0025) (−3 dB beamwidth),  φ1 = (φ0/2) √10.1875
/// C = 21 − 25 log φ1 − (Gmax − 17)
/// ```
///
/// RR Appendix 30 (Annex 5, Fig. 7bis) uses this pattern for the 60 cm
/// reference antenna of the Regions 1 and 3 Plan.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuBo1213 {
    peak_gain_dbi: f64,
    d_over_lambda: f64,
    component: PolarComponent,
    phi_m_deg: f64,
    g1_dbi: f64,
    phi_r_deg: f64,
    phi_0_deg: f64,
    phi_1_deg: f64,
    c_db: f64,
}

impl ItuBo1213 {
    /// Pattern for `D/λ` (equivalent diameter) and antenna efficiency.
    ///
    /// ```
    /// use antenna_r::{GainPattern, parabolic::itu::{ItuBo1213, PolarComponent}};
    /// // 60 cm WRC-03 reference antenna: D/λ = 23.4, η = 0.65.
    /// let p = ItuBo1213::new(23.4, 0.65, PolarComponent::CoPolar).unwrap();
    /// assert!((p.peak_gain_dbi() - 35.5).abs() < 0.05);
    /// assert!((p.phi_m_deg() - 3.98).abs() < 0.01);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for `D/λ < 11` (or non-finite) or an efficiency outside `(0, 1]`.
    pub fn new(d_over_lambda: f64, efficiency: f64, component: PolarComponent) -> Result<Self> {
        let dl = require_d_over_lambda(d_over_lambda)?;
        if !(dl >= BO1213_MIN_D_OVER_LAMBDA && dl.is_finite()) {
            return Err(invalid(
                "d_over_lambda",
                format!("BO.1213 is valid for finite D/λ ≥ 11, got {dl}"),
            ));
        }
        let eta = require_efficiency("efficiency", efficiency)?;
        let gmax = rf_r::db::to_db(eta * (std::f64::consts::PI * dl).powi(2));
        let phi_r_deg = 95.0 / dl;
        let g1_dbi = 29.0 - 25.0 * phi_r_deg.log10();
        let phi_0_deg = 2.0 / dl * (3.0f64 / 0.0025).sqrt();
        let phi_1_deg = phi_0_deg / 2.0 * 10.1875f64.sqrt();
        Ok(Self {
            peak_gain_dbi: gmax,
            d_over_lambda: dl,
            component,
            phi_m_deg: phi_m_deg(dl, gmax, g1_dbi),
            g1_dbi,
            phi_r_deg,
            phi_0_deg,
            phi_1_deg,
            c_db: 21.0 - 25.0 * phi_1_deg.log10() - (gmax - 17.0),
        })
    }

    /// Pattern of a dish of diameter `diameter_m` at `freq_hz`.
    ///
    /// # Errors
    /// As [`ItuBo1213::new`], plus a non-positive diameter or frequency.
    pub fn from_dish(diameter_m: f64, freq_hz: f64, efficiency: f64, component: PolarComponent) -> Result<Self> {
        Self::new(dish_d_over_lambda(diameter_m, freq_hz)?, efficiency, component)
    }

    /// `φm` (degrees), end of the co-polar main lobe.
    #[must_use]
    pub fn phi_m_deg(&self) -> f64 {
        self.phi_m_deg
    }

    /// `φ0` (degrees), the −3 dB beamwidth used by the cross-polar pattern.
    #[must_use]
    pub fn beamwidth_deg(&self) -> f64 {
        self.phi_0_deg
    }

    /// Constant `C` (dB) of the cross-polar pattern (negative, as the
    /// recommendation requires).
    #[must_use]
    pub fn cross_polar_c_db(&self) -> f64 {
        self.c_db
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        if phi_deg.is_nan() {
            return f64::NAN;
        }
        match self.component {
            PolarComponent::CoPolar => self.co_polar_dbi(phi_deg),
            PolarComponent::CrossPolar => self.cross_polar_dbi(phi_deg),
        }
    }

    fn co_polar_dbi(&self, phi_deg: f64) -> f64 {
        if phi_deg < self.phi_m_deg {
            main_lobe_parabola_dbi(self.peak_gain_dbi, self.d_over_lambda, phi_deg)
        } else if phi_deg < self.phi_r_deg {
            self.g1_dbi
        } else if phi_deg < BO1213_PHI_B_DEG {
            29.0 - 25.0 * phi_deg.log10()
        } else {
            back_levels_dbi(phi_deg)
        }
    }

    fn cross_polar_dbi(&self, phi_deg: f64) -> f64 {
        let g = self.peak_gain_dbi;
        let phi0 = self.phi_0_deg;
        if phi_deg < 0.25 * phi0 {
            g - 25.0
        } else if phi_deg < 0.44 * phi0 {
            g - 25.0 + 8.0 * (phi_deg - 0.25 * phi0) / (0.19 * phi0)
        } else if phi_deg < phi0 {
            g - 17.0
        } else if phi_deg < self.phi_1_deg {
            g - 17.0 + self.c_db * (phi_deg - phi0) / (self.phi_1_deg - phi0)
        } else if phi_deg < BO1213_PHI_2_DEG {
            21.0 - 25.0 * phi_deg.log10()
        } else {
            back_levels_dbi(phi_deg)
        }
    }
}

/// −5 dBi up to 70°, 0 dBi beyond (both BO.1213 components).
fn back_levels_dbi(phi_deg: f64) -> f64 {
    if phi_deg < BO1213_BACK_DEG { -5.0 } else { 0.0 }
}

impl GainPattern for ItuBo1213 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    /// Maximum of the selected component: `Gmax` (co-polar) or
    /// `Gmax − 17` (cross-polar).
    fn peak_gain_dbi(&self) -> f64 {
        match self.component {
            PolarComponent::CoPolar => self.peak_gain_dbi,
            PolarComponent::CrossPolar => self.peak_gain_dbi - 17.0,
        }
    }
}

/// On-axis gain taken by BO.1295 and RR Appendix 30A for the 5 m (17 GHz)
/// and 6 m (14 GHz) planning antennas (dBi).
pub const BO1295_PEAK_GAIN_DBI: f64 = 57.0;

/// ITU-R BO.1295-0 reference off-axis e.i.r.p. of BSS feeder-link
/// transmitting earth stations (14 and 17 GHz, Regions 1 and 3), as a gain
/// pattern: the off-axis e.i.r.p. relative to the on-axis e.i.r.p. `E` is
/// added to the peak gain. θ in degrees:
///
/// ```text
/// Co-polar (Curve A′):                 Cross-polar (Curve B′):
///   0                    θ ≤ 0.1°        −35                 θ ≤ 1.91°
///   −21 − 20 log θ       0.1° < θ ≤ 0.32°  −28 − 25 log θ    1.91° < θ ≤ 36.31°
///   −5.7 − 53.2 θ²       0.32° < θ ≤ 0.54° −67               θ > 36.31°
///   −28 − 25 log θ       0.54° < θ ≤ 36.31°
///   −67                  θ > 36.31°
/// ```
///
/// These are the WRC-97 curves of RR Appendix 30A (Annex 3, Fig. A); the
/// original WARC Orb-88 curves are in [`Ap30aEarthTransmit`](super::Ap30aEarthTransmit).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ItuBo1295 {
    peak_gain_dbi: f64,
    component: PolarComponent,
}

impl ItuBo1295 {
    /// Pattern with on-axis gain `peak_gain_dbi` (the recommendation takes
    /// [`BO1295_PEAK_GAIN_DBI`]).
    ///
    /// ```
    /// use antenna_r::parabolic::itu::{BO1295_PEAK_GAIN_DBI, ItuBo1295, PolarComponent};
    /// let p = ItuBo1295::new(BO1295_PEAK_GAIN_DBI, PolarComponent::CoPolar).unwrap();
    /// assert_eq!(p.gain_at_deg(50.0), 57.0 - 67.0);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain.
    pub fn new(peak_gain_dbi: f64, component: PolarComponent) -> Result<Self> {
        Ok(Self {
            peak_gain_dbi: require_finite("peak_gain_dbi", peak_gain_dbi)?,
            component,
        })
    }

    /// Gain at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, theta_deg: f64) -> f64 {
        self.peak_gain_dbi + bo1295_relative_db(self.component, theta_deg)
    }
}

/// BO.1295-0 off-axis e.i.r.p. relative to the on-axis value (dB).
pub(super) fn bo1295_relative_db(component: PolarComponent, theta_deg: f64) -> f64 {
    match component {
        PolarComponent::CoPolar => {
            if theta_deg <= 0.1 {
                0.0
            } else if theta_deg <= 0.32 {
                -21.0 - 20.0 * theta_deg.log10()
            } else if theta_deg <= 0.54 {
                -5.7 - 53.2 * theta_deg * theta_deg
            } else {
                earth_station_tail_db(theta_deg)
            }
        }
        PolarComponent::CrossPolar => {
            if theta_deg <= 1.91 {
                -35.0
            } else {
                earth_station_tail_db(theta_deg)
            }
        }
    }
}

/// `−28 − 25 log θ` up to 36.31°, then −67 dB (BO.1295 both components).
fn earth_station_tail_db(theta_deg: f64) -> f64 {
    if theta_deg <= 36.31 {
        -28.0 - 25.0 * theta_deg.log10()
    } else if theta_deg.is_nan() {
        f64::NAN
    } else {
        -67.0
    }
}

impl GainPattern for ItuBo1295 {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    /// Maximum of the selected component: the on-axis gain (co-polar) or
    /// 35 dB below it (cross-polar).
    fn peak_gain_dbi(&self) -> f64 {
        match self.component {
            PolarComponent::CoPolar => self.peak_gain_dbi,
            PolarComponent::CrossPolar => self.peak_gain_dbi - 35.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    #[test]
    fn bo1213_constants() {
        assert_close(BO1213_PHI_B_DEG, 10f64.powf(34.0 / 25.0), 1e-12);
        assert_close(BO1213_PHI_2_DEG, 10f64.powf(26.0 / 25.0), 1e-12);
    }

    #[test]
    fn bo1213_matches_the_60_cm_example() {
        let co = ItuBo1213::new(23.4, 0.65, PolarComponent::CoPolar).unwrap();
        assert_close(co.peak_gain_dbi(), 35.5, 0.05);
        assert_close(co.phi_m_deg(), 3.98, 0.01);
        assert_close(co.phi_r_deg, 4.06, 0.01);
        assert_close(co.g1_dbi, 13.78, 0.01);
        assert_close(co.beamwidth_deg(), 2.96, 0.01);
        assert_close(co.phi_1_deg, 4.73, 0.01);
        assert_close(co.cross_polar_c_db(), -14.36, 0.05);
        let cross = ItuBo1213::new(23.4, 0.65, PolarComponent::CrossPolar).unwrap();
        assert_close(cross.peak_gain_dbi(), co.peak_gain_dbi() - 17.0, 1e-12);
    }

    #[test]
    fn bo1213_matches_the_45_cm_example() {
        let p = ItuBo1213::new(18.3, 0.65, PolarComponent::CoPolar).unwrap();
        assert_close(p.peak_gain_dbi(), 33.3, 0.05);
        assert_close(p.phi_m_deg(), 5.15, 0.01);
        assert_close(p.g1_dbi, 11.12, 0.01);
        assert_close(p.beamwidth_deg(), 3.79, 0.01);
        assert_close(p.phi_1_deg, 6.04, 0.01);
        assert_close(p.cross_polar_c_db(), -14.83, 0.05);
    }

    #[test]
    fn bo1213_co_polar_segments() {
        let p = ItuBo1213::new(23.4, 0.65, PolarComponent::CoPolar).unwrap();
        let g = p.peak_gain_dbi();
        assert_eq!(p.gain_dbi(Direction::BORESIGHT), g);
        assert_close(p.gain_at_deg(2.0), g - 2.5e-3 * (23.4f64 * 2.0).powi(2), 1e-12);
        assert_eq!(p.gain_at_deg(4.0), p.g1_dbi);
        assert_close(p.gain_at_deg(10.0), 4.0, 1e-12);
        assert_eq!(p.gain_at_deg(23.0), -5.0);
        assert_close(p.gain_at_deg(22.9), 29.0 - 25.0 * 22.9f64.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(70.0), 0.0);
        assert_eq!(p.gain_at_deg(180.0), 0.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn bo1213_cross_polar_segments() {
        let p = ItuBo1213::new(23.4, 0.65, PolarComponent::CrossPolar).unwrap();
        let g = p.peak_gain_dbi;
        let phi0 = p.beamwidth_deg();
        assert_eq!(p.gain_at_deg(0.0), g - 25.0);
        assert_close(p.gain_at_deg(0.25 * phi0), g - 25.0, 1e-12);
        assert_close(p.gain_at_deg(0.44 * phi0 - 1e-12), g - 17.0, 1e-9);
        assert_eq!(p.gain_at_deg(0.5 * phi0), g - 17.0);
        assert_close(p.gain_at_deg(phi0), g - 17.0, 1e-12);
        // Linear ramp reaches 21 − 25 log φ1 at φ1.
        let phi1 = p.phi_1_deg;
        assert_close(p.gain_at_deg(phi1 - 1e-12), 21.0 - 25.0 * phi1.log10(), 1e-9);
        assert_close(p.gain_at_deg(8.0), 21.0 - 25.0 * 8f64.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(11.0), -5.0);
        assert_eq!(p.gain_at_deg(90.0), 0.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn bo1213_validation_and_from_dish() {
        for bad in [10.9, 0.0, f64::INFINITY] {
            assert!(matches!(
                ItuBo1213::new(bad, 0.65, PolarComponent::CoPolar),
                Err(AntennaError::InvalidParameter {
                    name: "d_over_lambda",
                    ..
                })
            ));
        }
        assert!(matches!(
            ItuBo1213::new(20.0, 0.0, PolarComponent::CoPolar),
            Err(AntennaError::InvalidParameter { name: "efficiency", .. })
        ));
        let p = ItuBo1213::from_dish(0.6, 11.7e9, 0.65, PolarComponent::CoPolar).unwrap();
        assert_close(p.d_over_lambda, 0.6 / rf_r::fspl::wavelength_m(11.7e9), 1e-12);
        assert!(ItuBo1213::from_dish(0.6, 0.0, 0.65, PolarComponent::CoPolar).is_err());
    }

    #[test]
    fn bo1295_co_polar_curve() {
        let p = ItuBo1295::new(BO1295_PEAK_GAIN_DBI, PolarComponent::CoPolar).unwrap();
        assert_eq!(p.gain_dbi(Direction::BORESIGHT), 57.0);
        assert_eq!(p.gain_at_deg(0.1), 57.0);
        assert_close(p.gain_at_deg(0.2), 57.0 - 21.0 - 20.0 * 0.2f64.log10(), 1e-12);
        assert_close(p.gain_at_deg(0.5), 57.0 - 5.7 - 53.2 * 0.25, 1e-12);
        assert_close(p.gain_at_deg(1.0), 57.0 - 28.0, 1e-12);
        assert_close(p.gain_at_deg(36.31), 57.0 - 28.0 - 25.0 * 36.31f64.log10(), 1e-12);
        assert_eq!(p.gain_at_deg(36.32), 57.0 - 67.0);
        assert_eq!(p.peak_gain_dbi(), 57.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn bo1295_cross_polar_curve_and_validation() {
        let p = ItuBo1295::new(57.0, PolarComponent::CrossPolar).unwrap();
        assert_eq!(p.gain_at_deg(1.0), 22.0);
        assert_close(p.gain_at_deg(10.0), 57.0 - 53.0, 1e-12);
        assert_eq!(p.gain_at_deg(40.0), -10.0);
        assert_eq!(p.peak_gain_dbi(), 22.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
        assert!(ItuBo1295::new(f64::NAN, PolarComponent::CoPolar).is_err());
    }
}
