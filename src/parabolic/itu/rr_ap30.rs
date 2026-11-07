//! Reference patterns of the Radio Regulations Appendices 30 (BSS downlink
//! Plans, Annex 5 § 3.7.2 and § 3.13.3) and 30A (feeder-link Plans,
//! Annex 3 § 3.5.3 and § 3.7.3), Edition 2024, often called APS30 and
//! APS30A patterns.
//!
//! Relative angles are `x = φ/φ0` with `φ0` the full half-power beamwidth
//! in the plane of interest (degrees). Relative gains are in dB below the
//! on-axis gain; "Curve C" (minus the on-axis gain, i.e. 0 dBi) bounds the
//! satellite patterns and the community-reception earth-station pattern
//! from below.

use super::{BeamEllipse, PolarComponent, bss::bo1295_relative_db, off_axis_deg, plane_rad};
use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite};
use crate::pattern::GainPattern;

/// Largest full beamwidth accepted (degrees).
const MAX_BEAMWIDTH_DEG: f64 = 180.0;
/// Constant of `Gm = 27 843 / (a b)` (a, b in degrees, 55 % efficiency).
const AP30_GAIN_BEAMWIDTH_PRODUCT: f64 = 27_843.0;

/// Earth-station receiving pattern of RR Appendix 30 (Annex 5 § 3.7.2).
/// The WRC-97 Regions 1 and 3 pattern (Fig. 7bis) is [`ItuBo1213`](super::ItuBo1213).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ap30EarthKind {
    /// Fig. 7, Curves A (co-polar) and B (cross-polar): individual reception
    /// in the original 1977 Regions 1 and 3 Plan.
    Regions1And3Individual,
    /// Fig. 7, Curve A′ then C (co-polar) and B (cross-polar): community
    /// reception in the original 1977 Regions 1 and 3 Plan.
    Regions1And3Community,
    /// Fig. 8, Curves A and B: Region 2 (`φ0` = 1.7° for planning).
    Region2,
}

/// RR Appendix 30 reference receiving earth-station pattern (Annex 5,
/// Figs. 7 and 8), `x = φ/φ0`, relative gain in dB:
///
/// ```text
/// Fig. 7 Curve A (individual):   0 (x ≤ 0.25); −12 x² (≤ 0.707);
///   −(9.0 + 20 log x) (≤ 1.26); −(8.5 + 25 log x) (≤ 9.55); −33 beyond
/// Fig. 7 Curve A′ (community):   0 (x ≤ 0.25); −12 x² (≤ 0.86);
///   −(10.5 + 25 log x) beyond, down to Curve C (−Gon)
/// Fig. 7 Curve B (cross):        −25 (x ≤ 0.25); −(30 + 40 log|x − 1|) (≤ 0.44);
///   −20 (≤ 1.4); −(30 + 25 log|x − 1|) (≤ 2); −30 until it meets the
///   co-polar curve, then the co-polar curve
/// Fig. 8 Curve A (Region 2):     0 (x ≤ 0.25); −12 x² (≤ 1.13);
///   −(14 + 25 log x) (≤ 14.7); −43.2 (≤ 35); −(85.2 − 27.2 log x) (≤ 45.1);
///   −40.2 (≤ 70); −(−55.2 + 51.7 log x) (≤ 80); −43.2 to 180°
/// Fig. 8 Curve B (Region 2):     −25 (x ≤ 0.25); −(30 + 40 log|x − 1|) (≤ 0.44);
///   −20 (≤ 1.28); −(17.3 + 25 log x) (≤ 3.22); −30 until it meets the
///   co-polar curve, then the co-polar curve
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ap30EarthReceive {
    peak_gain_dbi: f64,
    phi0_deg: f64,
    kind: Ap30EarthKind,
    component: PolarComponent,
}

impl Ap30EarthReceive {
    /// Pattern with on-axis gain `peak_gain_dbi` and full half-power
    /// beamwidth `beamwidth_rad` (`φ0`).
    ///
    /// ```
    /// use antenna_r::parabolic::itu::{Ap30EarthKind, Ap30EarthReceive, PolarComponent};
    /// let p = Ap30EarthReceive::new(37.0, 2f64.to_radians(), Ap30EarthKind::Regions1And3Individual,
    ///     PolarComponent::CoPolar).unwrap();
    /// assert_eq!(p.gain_at_deg(30.0), 37.0 - 33.0);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain or a beamwidth outside `(0, 180°)`.
    pub fn new(peak_gain_dbi: f64, beamwidth_rad: f64, kind: Ap30EarthKind, component: PolarComponent) -> Result<Self> {
        let peak = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let phi0_deg = beamwidth_rad.to_degrees();
        if !(phi0_deg > 0.0 && phi0_deg < MAX_BEAMWIDTH_DEG) {
            return Err(invalid(
                "beamwidth_rad",
                format!("must be in (0, π), got {beamwidth_rad}"),
            ));
        }
        Ok(Self {
            peak_gain_dbi: peak,
            phi0_deg,
            kind,
            component,
        })
    }

    /// Relative gain (dB, ≤ 0) at an off-axis angle in degrees.
    #[must_use]
    pub fn relative_gain_db(&self, phi_deg: f64) -> f64 {
        if phi_deg.is_nan() {
            return f64::NAN;
        }
        let x = phi_deg / self.phi0_deg;
        let co = match self.kind {
            Ap30EarthKind::Regions1And3Individual => fig7_individual_co_db(x),
            Ap30EarthKind::Regions1And3Community => fig7_community_co_db(x).max(-self.peak_gain_dbi),
            Ap30EarthKind::Region2 => fig8_co_db(x),
        };
        match self.component {
            PolarComponent::CoPolar => co,
            PolarComponent::CrossPolar => match self.kind {
                Ap30EarthKind::Region2 => fig8_cross_db(x, co),
                _ => fig7_cross_db(x, co),
            },
        }
    }

    /// Gain (dBi) at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64) -> f64 {
        self.peak_gain_dbi + self.relative_gain_db(phi_deg)
    }
}

/// Fig. 7 Curve A.
fn fig7_individual_co_db(x: f64) -> f64 {
    if x <= 0.25 {
        0.0
    } else if x <= 0.707 {
        -12.0 * x * x
    } else if x <= 1.26 {
        -(9.0 + 20.0 * x.log10())
    } else if x <= 9.55 {
        -(8.5 + 25.0 * x.log10())
    } else {
        -33.0
    }
}

/// Fig. 7 Curve A′ (before Curve C).
fn fig7_community_co_db(x: f64) -> f64 {
    if x <= 0.25 {
        0.0
    } else if x <= 0.86 {
        -12.0 * x * x
    } else {
        -(10.5 + 25.0 * x.log10())
    }
}

/// Fig. 7 Curve B; `co` is the co-polar relative gain at `x`.
fn fig7_cross_db(x: f64, co: f64) -> f64 {
    if x <= 0.25 {
        -25.0
    } else if x <= 0.44 {
        -(30.0 + 40.0 * (1.0 - x).log10())
    } else if x <= 1.4 {
        -20.0
    } else if x <= 2.0 {
        -(30.0 + 25.0 * (x - 1.0).log10())
    } else {
        co.min(-30.0)
    }
}

/// Fig. 8 Curve A.
fn fig8_co_db(x: f64) -> f64 {
    let log_x = x.log10();
    if x <= 0.25 {
        0.0
    } else if x <= 1.13 {
        -12.0 * x * x
    } else if x <= 14.7 {
        -(14.0 + 25.0 * log_x)
    } else if x <= 35.0 {
        -43.2
    } else if x <= 45.1 {
        -(85.2 - 27.2 * log_x)
    } else if x <= 70.0 {
        -40.2
    } else if x <= 80.0 {
        -(-55.2 + 51.7 * log_x)
    } else {
        -43.2
    }
}

/// Fig. 8 Curve B; `co` is the co-polar relative gain at `x`.
fn fig8_cross_db(x: f64, co: f64) -> f64 {
    if x <= 0.25 {
        -25.0
    } else if x <= 0.44 {
        -(30.0 + 40.0 * (1.0 - x).log10())
    } else if x <= 1.28 {
        -20.0
    } else if x <= 3.22 {
        -(17.3 + 25.0 * x.log10())
    } else {
        co.min(-30.0)
    }
}

impl GainPattern for Ap30EarthReceive {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    /// Maximum of the selected component: the on-axis gain (co-polar) or
    /// 20 dB below it (cross-polar).
    fn peak_gain_dbi(&self) -> f64 {
        match self.component {
            PolarComponent::CoPolar => self.peak_gain_dbi,
            PolarComponent::CrossPolar => self.peak_gain_dbi - 20.0,
        }
    }
}

/// Satellite antenna pattern of RR Appendix 30 (transmitting, Annex 5
/// § 3.13.3) or Appendix 30A (receiving, Annex 3 § 3.7.3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ap30SatelliteKind {
    /// AP30 Fig. 9: transmitting antennas, Regions 1 and 3.
    TransmitRegions1And3,
    /// AP30 Fig. 10: transmitting antennas, Region 2.
    TransmitRegion2,
    /// AP30 Fig. 11B (transmitting, Regions 1 and 3) and AP30A Fig. C
    /// (receiving): fast roll-off in the main beam, beamlet 0.6°.
    FastRollOffRegions1And3,
    /// AP30 Fig. 11A: fast roll-off transmitting antennas, Region 2,
    /// beamlet 0.8°.
    FastRollOffRegion2,
    /// AP30 Fig. 13 (ITU-R BO.1445): improved fast roll-off transmitting
    /// antennas, Regions 1 and 3, beamlet 0.6°.
    ImprovedFastRollOffRegions1And3,
    /// AP30A Fig. B, Curves A and B: receiving antennas, WARC Orb-88 Plan.
    ReceiveOrb88,
    /// AP30A Fig. B, Curves A′ and B′ (ITU-R BO.1296): receiving antennas,
    /// WRC-97 replanning.
    ReceiveWrc97,
}

impl Ap30SatelliteKind {
    /// Beamlet half-power beamwidth `Bmin` (degrees) of the fast roll-off
    /// patterns.
    fn beamlet_deg(self) -> Option<f64> {
        match self {
            Self::FastRollOffRegions1And3 | Self::ImprovedFastRollOffRegions1And3 => Some(0.6),
            Self::FastRollOffRegion2 => Some(0.8),
            _ => None,
        }
    }

    /// Peak of the cross-polar relative pattern (dB).
    fn cross_polar_peak_db(self) -> f64 {
        match self {
            Self::TransmitRegions1And3 | Self::ImprovedFastRollOffRegions1And3 => -33.0,
            Self::ReceiveWrc97 => -35.0,
            _ => -30.0,
        }
    }
}

/// RR Appendix 30 / 30A satellite antenna reference patterns for circular
/// or elliptical beams, `x = φ/φ0`, relative gain in dB, all bounded below
/// by Curve C (`−Gon`):
///
/// ```text
/// AP30 Fig. 9   co:    −12 x² (x ≤ 1.58); −30 (≤ 3.16); −(17.5 + 25 log x)
///               cross: −(40 + 40 log|x − 1|) (x ≤ 0.33); −33 (≤ 1.67);
///                      −(40 + 40 log|x − 1|) beyond
/// AP30 Fig. 10  co:    −12 x² (x ≤ 1.45); −(22 + 20 log x)
///               cross: −30 (x ≤ 2.51), then the co-polar curve
/// Fast roll-off (AP30 Figs. 11A/11B, AP30A Fig. C), beamlet Bmin:
///               co:    −12 x² (x ≤ 0.5); −12 ((x − x0)/(Bmin/φ0))²
///                      (≤ 1.45 Bmin/φ0 + x0); −25.23 (≤ 1.45); −(22 + 20 log x)
///                      x0 = 0.5 (1 − Bmin/φ0)
///               cross: −30 (x < 2.51), then the co-polar curve
/// AP30 Fig. 13  co:    min(fast roll-off with Bmin = 0.6°, Fig. 9 co)
///               cross: Fig. 9 cross
/// AP30A Fig. B  co:    −12 x² (x < 1.3); −17.5 − 25 log x
///   Orb-88      cross: −30 − 12 x² (x < 0.5); −33 (< 1.67); −40 − 40 log|x − 1|
///   WRC-97      cross: −35 (x < 1.75); −40 − 40 log|x − 1|
/// ```
///
/// `φ0` follows an ellipse with full beamwidths `a` (major axis, along
/// azimuth `orientation`) and `b`. The Plans relate peak gain and
/// beamwidths by `Gm = 27 843 / (a b)` (55 % efficiency), see
/// [`Ap30Satellite::peak_gain_from_beamwidths_dbi`]. The fast roll-off
/// patterns need `φ0 ≥ Bmin` (the minimum beamwidth of the Plans).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ap30Satellite {
    peak_gain_dbi: f64,
    beamwidth: BeamEllipse,
    kind: Ap30SatelliteKind,
    component: PolarComponent,
}

impl Ap30Satellite {
    /// Circular beam with full half-power beamwidth `beamwidth_rad`.
    ///
    /// ```
    /// use antenna_r::parabolic::itu::{Ap30Satellite, Ap30SatelliteKind, PolarComponent};
    /// let p = Ap30Satellite::circular(43.0, 1.2f64.to_radians(), Ap30SatelliteKind::TransmitRegions1And3,
    ///     PolarComponent::CoPolar).unwrap();
    /// assert!((p.gain_at_deg(0.6, 1.2) - 40.0).abs() < 1e-12); // −3 dB at φ0/2
    /// assert_eq!(p.gain_at_deg(3.0, 1.2), 13.0);                 // −30 dB plateau
    /// ```
    ///
    /// # Errors
    /// As [`Ap30Satellite::elliptical`].
    pub fn circular(
        peak_gain_dbi: f64,
        beamwidth_rad: f64,
        kind: Ap30SatelliteKind,
        component: PolarComponent,
    ) -> Result<Self> {
        Self::elliptical(peak_gain_dbi, beamwidth_rad, beamwidth_rad, 0.0, kind, component)
    }

    /// Elliptical beam with full half-power beamwidths along (`a`) and
    /// across (`b`) the major axis, which lies at azimuth `orientation_rad`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain or orientation, beamwidths outside `(0, 180°)`,
    /// a minor beamwidth larger than the major one, or (fast roll-off
    /// patterns) a minor beamwidth below the beamlet beamwidth.
    pub fn elliptical(
        peak_gain_dbi: f64,
        major_beamwidth_rad: f64,
        minor_beamwidth_rad: f64,
        orientation_rad: f64,
        kind: Ap30SatelliteKind,
        component: PolarComponent,
    ) -> Result<Self> {
        let peak = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        let beamwidth = BeamEllipse::new(
            ("major_beamwidth_rad", "minor_beamwidth_rad"),
            major_beamwidth_rad,
            minor_beamwidth_rad,
            orientation_rad,
            MAX_BEAMWIDTH_DEG,
        )?;
        if let Some(bmin) = kind.beamlet_deg()
            && beamwidth.minor_deg < bmin
        {
            return Err(invalid(
                "minor_beamwidth_rad",
                format!("fast roll-off patterns need φ0 ≥ {bmin}°, got {}°", beamwidth.minor_deg),
            ));
        }
        Ok(Self {
            peak_gain_dbi: peak,
            beamwidth,
            kind,
            component,
        })
    }

    /// `Gm = 10 log(27 843 / (a b))` dBi for full beamwidths `a`, `b`
    /// (AP30 Annex 5 § 3.13.1, AP30A Annex 3 § 3.7.1).
    ///
    /// ```
    /// use antenna_r::parabolic::itu::Ap30Satellite;
    /// let g = Ap30Satellite::peak_gain_from_beamwidths_dbi(1f64.to_radians(), 1f64.to_radians()).unwrap();
    /// assert!((g - 44.447).abs() < 1e-3);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for non-positive beamwidths.
    pub fn peak_gain_from_beamwidths_dbi(major_beamwidth_rad: f64, minor_beamwidth_rad: f64) -> Result<f64> {
        let a = crate::error::require_positive("major_beamwidth_rad", major_beamwidth_rad)?.to_degrees();
        let b = crate::error::require_positive("minor_beamwidth_rad", minor_beamwidth_rad)?.to_degrees();
        Ok(rf_r::db::to_db(AP30_GAIN_BEAMWIDTH_PRODUCT / (a * b)))
    }

    /// Full half-power beamwidth `φ0` (degrees) in the plane `φ = phi_rad`.
    #[must_use]
    pub fn beamwidth_deg(&self, phi_rad: f64) -> f64 {
        self.beamwidth.width_deg(phi_rad)
    }

    /// Relative gain (dB, ≤ 0) at off-axis angle `phi_deg` in a plane where
    /// the beamwidth is `phi0_deg`, Curve C included.
    #[must_use]
    pub fn relative_gain_db(&self, phi_deg: f64, phi0_deg: f64) -> f64 {
        if phi_deg.is_nan() {
            return f64::NAN;
        }
        let x = phi_deg / phi0_deg;
        let rel = match self.component {
            PolarComponent::CoPolar => self.co_polar_db(x, phi0_deg),
            PolarComponent::CrossPolar => self.cross_polar_db(x, phi0_deg),
        };
        rel.max(-self.peak_gain_dbi)
    }

    /// Gain (dBi) at off-axis angle `phi_deg` in a plane where the beamwidth
    /// is `phi0_deg`.
    #[must_use]
    pub fn gain_at_deg(&self, phi_deg: f64, phi0_deg: f64) -> f64 {
        self.peak_gain_dbi + self.relative_gain_db(phi_deg, phi0_deg)
    }

    fn co_polar_db(&self, x: f64, phi0_deg: f64) -> f64 {
        match self.kind {
            Ap30SatelliteKind::TransmitRegions1And3 => fig9_co_db(x),
            Ap30SatelliteKind::TransmitRegion2 => fig10_co_db(x),
            Ap30SatelliteKind::FastRollOffRegions1And3 => fast_roll_off_co_db(x, 0.6 / phi0_deg),
            Ap30SatelliteKind::FastRollOffRegion2 => fast_roll_off_co_db(x, 0.8 / phi0_deg),
            Ap30SatelliteKind::ImprovedFastRollOffRegions1And3 => {
                fast_roll_off_co_db(x, 0.6 / phi0_deg).min(fig9_co_db(x))
            }
            Ap30SatelliteKind::ReceiveOrb88 | Ap30SatelliteKind::ReceiveWrc97 => receive_co_db(x),
        }
    }

    fn cross_polar_db(&self, x: f64, phi0_deg: f64) -> f64 {
        match self.kind {
            Ap30SatelliteKind::TransmitRegions1And3 | Ap30SatelliteKind::ImprovedFastRollOffRegions1And3 => {
                fig9_cross_db(x)
            }
            Ap30SatelliteKind::ReceiveOrb88 => {
                if x < 0.5 {
                    -30.0 - 12.0 * x * x
                } else if x < 1.67 {
                    -33.0
                } else {
                    -40.0 - 40.0 * (x - 1.0).log10()
                }
            }
            Ap30SatelliteKind::ReceiveWrc97 => {
                if x < 1.75 {
                    -35.0
                } else {
                    -40.0 - 40.0 * (x - 1.0).log10()
                }
            }
            _ => {
                if x < 2.51 {
                    -30.0
                } else {
                    self.co_polar_db(x, phi0_deg)
                }
            }
        }
    }
}

/// AP30 Fig. 9 Curve A.
fn fig9_co_db(x: f64) -> f64 {
    if x <= 1.58 {
        -12.0 * x * x
    } else if x <= 3.16 {
        -30.0
    } else {
        -(17.5 + 25.0 * x.log10())
    }
}

/// AP30 Fig. 9 Curve B.
fn fig9_cross_db(x: f64) -> f64 {
    if x <= 0.33 || x > 1.67 {
        -(40.0 + 40.0 * (x - 1.0).abs().log10())
    } else {
        -33.0
    }
}

/// AP30 Fig. 10 Curve A.
fn fig10_co_db(x: f64) -> f64 {
    if x <= 1.45 {
        -12.0 * x * x
    } else {
        -(22.0 + 20.0 * x.log10())
    }
}

/// Fast roll-off co-polar curve (AP30 Figs. 11A, 11B and 13 `ΔG1`, AP30A
/// Fig. C); `beamlet_ratio = Bmin/φ0`.
fn fast_roll_off_co_db(x: f64, beamlet_ratio: f64) -> f64 {
    let x0 = 0.5 * (1.0 - beamlet_ratio);
    if x <= 0.5 {
        -12.0 * x * x
    } else if x <= 1.45f64.mul_add(beamlet_ratio, x0) {
        let r = (x - x0) / beamlet_ratio;
        -12.0 * r * r
    } else if x <= 1.45 {
        -25.23
    } else {
        -(22.0 + 20.0 * x.log10())
    }
}

/// AP30A Fig. B Curves A and A′.
fn receive_co_db(x: f64) -> f64 {
    if x < 1.3 {
        -12.0 * x * x
    } else {
        -17.5 - 25.0 * x.log10()
    }
}

impl GainPattern for Ap30Satellite {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let phi0 = self.beamwidth.width_deg(plane_rad(dir));
        self.gain_at_deg(off_axis_deg(dir), phi0)
    }
    /// Maximum of the selected component: the on-axis gain (co-polar) or
    /// the peak of the cross-polar curve.
    fn peak_gain_dbi(&self) -> f64 {
        match self.component {
            PolarComponent::CoPolar => self.peak_gain_dbi,
            PolarComponent::CrossPolar => self.peak_gain_dbi + self.kind.cross_polar_peak_db().max(-self.peak_gain_dbi),
        }
    }
}

/// Feeder-link Plan of RR Appendix 30A for transmitting earth stations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ap30aPlan {
    /// WARC Orb-88 Curves A and B (original 1988 Plan).
    Orb88,
    /// WRC-97 Curves A′ and B′, as in ITU-R BO.1295.
    Wrc97,
}

/// RR Appendix 30A (Annex 3 § 3.5.3, Fig. A) off-axis e.i.r.p. of
/// feeder-link transmitting earth stations, as a gain pattern (θ in degrees,
/// dB relative to the on-axis value):
///
/// ```text
/// Orb-88 co-polar (A):    0 (θ < 0.1°); −21 − 20 log θ (< 0.32°);
///                         −5.7 − 53.2 θ² (< 0.44°); −25 − 25 log θ (< 48°); −67
/// Orb-88 cross-polar (B): −30 (θ < 1.6°); −25 − 25 log θ (< 48°); −67
/// WRC-97 (A′, B′):        as ItuBo1295
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ap30aEarthTransmit {
    peak_gain_dbi: f64,
    plan: Ap30aPlan,
    component: PolarComponent,
}

impl Ap30aEarthTransmit {
    /// Pattern with on-axis gain `peak_gain_dbi` (the Plan takes
    /// [`BO1295_PEAK_GAIN_DBI`](super::BO1295_PEAK_GAIN_DBI)).
    ///
    /// ```
    /// use antenna_r::parabolic::itu::{Ap30aEarthTransmit, Ap30aPlan, PolarComponent};
    /// let p = Ap30aEarthTransmit::new(57.0, Ap30aPlan::Orb88, PolarComponent::CrossPolar).unwrap();
    /// assert_eq!(p.gain_at_deg(1.0), 27.0);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain.
    pub fn new(peak_gain_dbi: f64, plan: Ap30aPlan, component: PolarComponent) -> Result<Self> {
        Ok(Self {
            peak_gain_dbi: require_finite("peak_gain_dbi", peak_gain_dbi)?,
            plan,
            component,
        })
    }

    /// Relative off-axis e.i.r.p. (dB) at an off-axis angle in degrees.
    #[must_use]
    pub fn relative_gain_db(&self, theta_deg: f64) -> f64 {
        match self.plan {
            Ap30aPlan::Wrc97 => bo1295_relative_db(self.component, theta_deg),
            Ap30aPlan::Orb88 => orb88_earth_relative_db(self.component, theta_deg),
        }
    }

    /// Gain (dBi) at an off-axis angle in degrees.
    #[must_use]
    pub fn gain_at_deg(&self, theta_deg: f64) -> f64 {
        self.peak_gain_dbi + self.relative_gain_db(theta_deg)
    }
}

/// AP30A Fig. A Curves A and B (WARC Orb-88).
fn orb88_earth_relative_db(component: PolarComponent, theta_deg: f64) -> f64 {
    if theta_deg.is_nan() {
        return f64::NAN;
    }
    let side = |t: f64| -25.0 - 25.0 * t.log10();
    match component {
        PolarComponent::CoPolar => {
            if theta_deg < 0.1 {
                0.0
            } else if theta_deg < 0.32 {
                -21.0 - 20.0 * theta_deg.log10()
            } else if theta_deg < 0.44 {
                -5.7 - 53.2 * theta_deg * theta_deg
            } else if theta_deg < 48.0 {
                side(theta_deg)
            } else {
                -67.0
            }
        }
        PolarComponent::CrossPolar => {
            if theta_deg < 1.6 {
                -30.0
            } else if theta_deg < 48.0 {
                side(theta_deg)
            } else {
                -67.0
            }
        }
    }
}

impl GainPattern for Ap30aEarthTransmit {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        self.gain_at_deg(off_axis_deg(dir))
    }
    /// Maximum of the selected component: the on-axis gain (co-polar) or
    /// 30 dB (Orb-88) / 35 dB (WRC-97) below it (cross-polar).
    fn peak_gain_dbi(&self) -> f64 {
        match (self.component, self.plan) {
            (PolarComponent::CoPolar, _) => self.peak_gain_dbi,
            (PolarComponent::CrossPolar, Ap30aPlan::Orb88) => self.peak_gain_dbi - 30.0,
            (PolarComponent::CrossPolar, Ap30aPlan::Wrc97) => self.peak_gain_dbi - 35.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    fn earth(kind: Ap30EarthKind, component: PolarComponent) -> Ap30EarthReceive {
        Ap30EarthReceive::new(37.0, 2f64.to_radians(), kind, component).unwrap()
    }

    #[test]
    fn fig7_individual_co_polar_is_continuous() {
        let p = earth(Ap30EarthKind::Regions1And3Individual, PolarComponent::CoPolar);
        assert_eq!(p.relative_gain_db(0.5), 0.0);
        assert_close(p.relative_gain_db(1.0), -3.0, 1e-12);
        for (x, tol) in [(0.707, 0.01), (1.26, 0.02), (9.55, 0.02)] {
            let phi = 2.0 * x;
            assert_close(p.relative_gain_db(phi - 1e-9), p.relative_gain_db(phi + 1e-9), tol);
        }
        assert_close(p.relative_gain_db(2.0 * 3.0), -(8.5 + 25.0 * 3f64.log10()), 1e-12);
        assert_eq!(p.relative_gain_db(40.0), -33.0);
        assert_eq!(p.gain_dbi(Direction::BORESIGHT), 37.0);
        assert_eq!(p.peak_gain_dbi(), 37.0);
        assert!(p.gain_at_deg(f64::NAN).is_nan());
    }

    #[test]
    fn fig7_community_co_polar_stops_at_curve_c() {
        let p = earth(Ap30EarthKind::Regions1And3Community, PolarComponent::CoPolar);
        assert_eq!(p.relative_gain_db(0.4), 0.0);
        assert_close(p.relative_gain_db(1.0), -3.0, 1e-12);
        assert_close(p.relative_gain_db(4.0), -(10.5 + 25.0 * 2f64.log10()), 1e-12);
        // −(10.5 + 25 log x) reaches −37 dB at x ≈ 87: Curve C beyond.
        assert_eq!(p.relative_gain_db(179.0), -37.0);
        assert_eq!(p.gain_at_deg(179.0), 0.0);
    }

    #[test]
    fn fig7_cross_polar_segments() {
        let p = earth(Ap30EarthKind::Regions1And3Individual, PolarComponent::CrossPolar);
        assert_eq!(p.relative_gain_db(0.4), -25.0);
        assert_close(p.relative_gain_db(0.5 + 1e-12), -25.0, 0.01);
        assert_close(p.relative_gain_db(0.7), -(30.0 + 40.0 * 0.65f64.log10()), 1e-12);
        assert_eq!(p.relative_gain_db(2.0), -20.0);
        assert_close(p.relative_gain_db(3.6), -(30.0 + 25.0 * 0.8f64.log10()), 1e-12);
        assert_eq!(p.relative_gain_db(6.0), -30.0);
        // Meets the co-polar curve near x = 7.24, then follows it.
        assert_close(p.relative_gain_db(2.0 * 8.0), -(8.5 + 25.0 * 8f64.log10()), 1e-12);
        assert_eq!(p.relative_gain_db(30.0), -33.0);
        assert_eq!(p.peak_gain_dbi(), 17.0);
        let community = earth(Ap30EarthKind::Regions1And3Community, PolarComponent::CrossPolar);
        assert_eq!(community.relative_gain_db(6.0), -30.0);
    }

    #[test]
    fn fig8_region2_curves() {
        let co = earth(Ap30EarthKind::Region2, PolarComponent::CoPolar);
        let at = |x: f64| co.relative_gain_db(2.0 * x);
        assert_eq!(at(0.2), 0.0);
        assert_close(at(1.0), -12.0, 1e-12);
        for x in [1.13, 14.7, 35.0, 45.1, 70.0, 80.0] {
            assert_close(at(x - 1e-9), at(x + 1e-9), 0.05);
        }
        assert_close(at(5.0), -(14.0 + 25.0 * 5f64.log10()), 1e-12);
        assert_eq!(at(20.0), -43.2);
        assert_close(at(40.0), -(85.2 - 27.2 * 40f64.log10()), 1e-12);
        assert_eq!(at(60.0), -40.2);
        assert_close(at(75.0), -(-55.2 + 51.7 * 75f64.log10()), 1e-12);
        assert_eq!(at(85.0), -43.2);
        let cross = earth(Ap30EarthKind::Region2, PolarComponent::CrossPolar);
        let xat = |x: f64| cross.relative_gain_db(2.0 * x);
        assert_eq!(xat(0.1), -25.0);
        assert_close(xat(0.3), -(30.0 + 40.0 * 0.7f64.log10()), 1e-12);
        assert_eq!(xat(1.0), -20.0);
        assert_close(xat(2.0), -(17.3 + 25.0 * 2f64.log10()), 1e-12);
        assert_close(xat(1.28 + 1e-9), -20.0, 0.03);
        assert_close(xat(3.22 + 1e-9), -30.0, 0.01);
        assert_eq!(xat(3.5), -30.0);
        // The co-polar curve falls below −30 dB near x = 3.98.
        assert_close(xat(5.0), -(14.0 + 25.0 * 5f64.log10()), 1e-12);
        assert_close(xat(10.0), -(14.0 + 25.0 * 10f64.log10()), 1e-12);
    }

    #[test]
    fn earth_receive_validation() {
        for bw in [0.0, -1.0, std::f64::consts::PI, f64::NAN] {
            assert!(matches!(
                Ap30EarthReceive::new(37.0, bw, Ap30EarthKind::Region2, PolarComponent::CoPolar),
                Err(AntennaError::InvalidParameter {
                    name: "beamwidth_rad",
                    ..
                })
            ));
        }
        assert!(Ap30EarthReceive::new(f64::NAN, 0.03, Ap30EarthKind::Region2, PolarComponent::CoPolar).is_err());
    }

    fn sat(kind: Ap30SatelliteKind, component: PolarComponent) -> Ap30Satellite {
        Ap30Satellite::circular(43.0, 1.2f64.to_radians(), kind, component).unwrap()
    }

    #[test]
    fn fig9_transmit_regions_1_and_3() {
        let co = sat(Ap30SatelliteKind::TransmitRegions1And3, PolarComponent::CoPolar);
        let at = |x: f64| co.relative_gain_db(1.2 * x, 1.2);
        assert_close(at(0.5), -3.0, 1e-12);
        assert_close(at(1.58), -12.0 * 1.58 * 1.58, 1e-12);
        assert_eq!(at(2.0), -30.0);
        assert_close(at(3.16 + 1e-9), -30.0, 0.01);
        assert_close(at(10.0), -42.5, 1e-12);
        // Curve C: never below −Gon.
        assert_eq!(at(100.0), -43.0);
        assert_eq!(co.gain_dbi(Direction::BORESIGHT), 43.0);
        assert_eq!(co.peak_gain_dbi(), 43.0);
        let cross = sat(Ap30SatelliteKind::TransmitRegions1And3, PolarComponent::CrossPolar);
        let xat = |x: f64| cross.relative_gain_db(1.2 * x, 1.2);
        assert_eq!(xat(0.0), -40.0);
        assert_close(xat(0.33), -33.0, 0.05);
        assert_eq!(xat(1.0), -33.0);
        assert_close(xat(1.67 + 1e-9), -33.0, 0.05);
        assert_close(xat(2.0), -40.0, 1e-12);
        // Curve C: −Gon = −43 dB.
        assert_eq!(xat(3.0), -43.0);
        assert_eq!(cross.peak_gain_dbi(), 10.0);
        assert!(co.gain_at_deg(f64::NAN, 1.2).is_nan());
    }

    #[test]
    fn fig10_transmit_region_2() {
        let co = sat(Ap30SatelliteKind::TransmitRegion2, PolarComponent::CoPolar);
        let at = |x: f64| co.relative_gain_db(1.2 * x, 1.2);
        assert_close(at(1.45), -25.23, 1e-12);
        assert_close(at(1.45 + 1e-9), -25.23, 0.01);
        assert_close(at(5.0), -(22.0 + 20.0 * 5f64.log10()), 1e-12);
        let cross = sat(Ap30SatelliteKind::TransmitRegion2, PolarComponent::CrossPolar);
        assert_eq!(cross.relative_gain_db(1.2 * 2.0, 1.2), -30.0);
        assert_close(cross.relative_gain_db(1.2 * 5.0, 1.2), at(5.0), 1e-12);
        assert_eq!(cross.peak_gain_dbi(), 13.0);
    }

    #[test]
    fn fast_roll_off_is_continuous_and_steeper() {
        for (kind, bmin) in [
            (Ap30SatelliteKind::FastRollOffRegions1And3, 0.6),
            (Ap30SatelliteKind::FastRollOffRegion2, 0.8),
        ] {
            let phi0 = 2.4_f64;
            let p = Ap30Satellite::circular(46.0, phi0.to_radians(), kind, PolarComponent::CoPolar).unwrap();
            let at = |x: f64| p.relative_gain_db(phi0 * x, phi0);
            let ratio = bmin / phi0;
            let x0 = 0.5 * (1.0 - ratio);
            let end = 1.45 * ratio + x0;
            assert_close(at(0.5), -3.0, 1e-12);
            assert_close(at(0.5 + 1e-9), -3.0, 1e-6);
            assert_close(at(end), -25.23, 1e-9);
            assert_eq!(at(f64::midpoint(end, 1.45)), -25.23);
            assert_close(at(3.0), -(22.0 + 20.0 * 3f64.log10()), 1e-12);
            // Steeper than the Fig. 10 main lobe beyond x = 0.5.
            assert!(at(0.9) < -12.0 * 0.81);
            let cross = Ap30Satellite::circular(46.0, phi0.to_radians(), kind, PolarComponent::CrossPolar).unwrap();
            assert_eq!(cross.relative_gain_db(phi0, phi0), -30.0);
            assert_close(cross.relative_gain_db(phi0 * 3.0, phi0), at(3.0), 1e-12);
        }
    }

    #[test]
    fn improved_fast_roll_off_takes_the_lower_curve() {
        let phi0 = 1.2;
        let p = sat(
            Ap30SatelliteKind::ImprovedFastRollOffRegions1And3,
            PolarComponent::CoPolar,
        );
        let fast = |x: f64| fast_roll_off_co_db(x, 0.6 / phi0);
        for x in [0.3, 0.8, 1.2, 2.0, 3.0, 5.0] {
            assert_close(p.relative_gain_db(phi0 * x, phi0), fast(x).min(fig9_co_db(x)), 1e-12);
        }
        // Between 1.45 and 3.16 the Fig. 9 plateau (−30) is lower than the
        // fast roll-off side lobes.
        assert_eq!(p.relative_gain_db(phi0 * 2.0, phi0), -30.0);
        let cross = sat(
            Ap30SatelliteKind::ImprovedFastRollOffRegions1And3,
            PolarComponent::CrossPolar,
        );
        assert_eq!(cross.relative_gain_db(phi0, phi0), -33.0);
    }

    #[test]
    fn ap30a_receive_patterns() {
        let co = sat(Ap30SatelliteKind::ReceiveOrb88, PolarComponent::CoPolar);
        let at = |x: f64| co.relative_gain_db(1.2 * x, 1.2);
        assert_close(at(1.0), -12.0, 1e-12);
        assert_close(at(2.0), -17.5 - 25.0 * 2f64.log10(), 1e-12);
        let wrc_co = sat(Ap30SatelliteKind::ReceiveWrc97, PolarComponent::CoPolar);
        assert_eq!(wrc_co.relative_gain_db(2.4, 1.2), at(2.0));
        let orb = sat(Ap30SatelliteKind::ReceiveOrb88, PolarComponent::CrossPolar);
        let oat = |x: f64| orb.relative_gain_db(1.2 * x, 1.2);
        assert_eq!(oat(0.0), -30.0);
        assert_close(oat(0.4), -30.0 - 12.0 * 0.16, 1e-12);
        assert_eq!(oat(1.0), -33.0);
        assert_close(oat(2.0), -40.0, 1e-12);
        assert_close(oat(1.9), -40.0 - 40.0 * 0.9f64.log10(), 1e-12);
        assert_eq!(orb.peak_gain_dbi(), 13.0);
        let wrc = sat(Ap30SatelliteKind::ReceiveWrc97, PolarComponent::CrossPolar);
        let wat = |x: f64| wrc.relative_gain_db(1.2 * x, 1.2);
        assert_eq!(wat(1.0), -35.0);
        assert_close(wat(1.75), -35.0, 0.01);
        assert_close(wat(2.0), -40.0, 1e-12);
        assert_eq!(wrc.peak_gain_dbi(), 8.0);
    }

    #[test]
    fn satellite_elliptical_beam_and_peak_gain() {
        let deg = 1f64.to_radians();
        let p = Ap30Satellite::elliptical(
            40.0,
            2.0 * deg,
            deg,
            0.0,
            Ap30SatelliteKind::TransmitRegion2,
            PolarComponent::CoPolar,
        )
        .unwrap();
        assert_close(p.beamwidth_deg(0.0), 2.0, 1e-12);
        assert_close(p.beamwidth_deg(std::f64::consts::FRAC_PI_2), 1.0, 1e-12);
        let along_major = p.gain_dbi(Direction::from_theta_phi(deg, 0.0));
        let along_minor = p.gain_dbi(Direction::from_theta_phi(0.5 * deg, std::f64::consts::FRAC_PI_2));
        assert_close(along_major, 37.0, 1e-9);
        assert_close(along_minor, 37.0, 1e-9);
        let g = Ap30Satellite::peak_gain_from_beamwidths_dbi(2.0 * deg, deg).unwrap();
        assert_close(g, 10.0 * (27_843.0f64 / 2.0).log10(), 1e-12);
        assert!(Ap30Satellite::peak_gain_from_beamwidths_dbi(0.0, deg).is_err());
        assert!(Ap30Satellite::peak_gain_from_beamwidths_dbi(deg, -deg).is_err());
        // Low on-axis gain: cross-polar peak limited by Curve C.
        let low = Ap30Satellite::circular(
            20.0,
            deg,
            Ap30SatelliteKind::TransmitRegion2,
            PolarComponent::CrossPolar,
        )
        .unwrap();
        assert_eq!(low.peak_gain_dbi(), 0.0);
    }

    #[test]
    fn satellite_validation() {
        let deg = 1f64.to_radians();
        let kind = Ap30SatelliteKind::TransmitRegion2;
        assert!(Ap30Satellite::circular(f64::NAN, deg, kind, PolarComponent::CoPolar).is_err());
        assert!(matches!(
            Ap30Satellite::circular(40.0, 0.0, kind, PolarComponent::CoPolar),
            Err(AntennaError::InvalidParameter {
                name: "major_beamwidth_rad",
                ..
            })
        ));
        assert!(matches!(
            Ap30Satellite::circular(
                40.0,
                0.5 * deg,
                Ap30SatelliteKind::FastRollOffRegions1And3,
                PolarComponent::CoPolar
            ),
            Err(AntennaError::InvalidParameter {
                name: "minor_beamwidth_rad",
                ..
            })
        ));
        assert!(
            Ap30Satellite::circular(
                40.0,
                0.7 * deg,
                Ap30SatelliteKind::FastRollOffRegion2,
                PolarComponent::CoPolar
            )
            .is_err()
        );
    }

    #[test]
    fn ap30a_earth_transmit_orb88_and_wrc97() {
        let co = Ap30aEarthTransmit::new(57.0, Ap30aPlan::Orb88, PolarComponent::CoPolar).unwrap();
        assert_eq!(co.relative_gain_db(0.05), 0.0);
        assert_close(co.relative_gain_db(0.2), -21.0 - 20.0 * 0.2f64.log10(), 1e-12);
        assert_close(co.relative_gain_db(0.4), -5.7 - 53.2 * 0.16, 1e-12);
        assert_close(co.relative_gain_db(10.0), -50.0, 1e-12);
        assert_eq!(co.relative_gain_db(48.0), -67.0);
        assert_eq!(co.gain_dbi(Direction::BORESIGHT), 57.0);
        assert_eq!(co.peak_gain_dbi(), 57.0);
        assert!(co.gain_at_deg(f64::NAN).is_nan());
        let cross = Ap30aEarthTransmit::new(57.0, Ap30aPlan::Orb88, PolarComponent::CrossPolar).unwrap();
        assert_eq!(cross.relative_gain_db(1.0), -30.0);
        assert_close(cross.relative_gain_db(10.0), -50.0, 1e-12);
        assert_eq!(cross.relative_gain_db(60.0), -67.0);
        assert_eq!(cross.peak_gain_dbi(), 27.0);
        let wrc = Ap30aEarthTransmit::new(57.0, Ap30aPlan::Wrc97, PolarComponent::CrossPolar).unwrap();
        assert_eq!(wrc.relative_gain_db(1.0), -35.0);
        assert_eq!(wrc.peak_gain_dbi(), 22.0);
        let wrc_co = Ap30aEarthTransmit::new(57.0, Ap30aPlan::Wrc97, PolarComponent::CoPolar).unwrap();
        assert_close(wrc_co.relative_gain_db(1.0), -28.0, 1e-12);
        assert!(Ap30aEarthTransmit::new(f64::INFINITY, Ap30aPlan::Wrc97, PolarComponent::CoPolar).is_err());
    }
}
