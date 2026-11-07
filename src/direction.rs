//! Directions in the antenna frame and the angle conventions used by pattern
//! files and positioners.
//!
//! Antenna frame: right-handed, boresight `+Z`, reference axis `+X`.
//! `θ` is the angle from `+Z` in `[0, π]`; `φ` is the azimuth in the X-Y
//! plane from `+X` towards `+Y` in `[0, 2π)`; `u = sinθ cosφ`,
//! `v = sinθ sinφ`, `w = cosθ`. These are the conventions of
//! [`geometry_r::satellite::ScanAngles`] and of TICRA GRASP.

use geometry_r::math::{Vec3, wrap_two_pi};
use geometry_r::satellite::ScanAngles;
use std::f64::consts::FRAC_PI_2;

use crate::error::{Result, invalid};

/// Tolerance on `u² + v² − 1` when building a direction from direction cosines.
const UV_TOLERANCE: f64 = 1e-12;

/// A unit direction in the antenna frame. `Copy`, 24 bytes.
///
/// Constructors that take user input validate it; the angle-based
/// constructors accept any value and propagate `NaN` for non-finite input.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction {
    v: Vec3,
}

impl Direction {
    /// The boresight `+Z`.
    pub const BORESIGHT: Self = Self { v: Vec3::Z };

    /// Wraps a vector that is already unit length (internal fast path).
    #[inline]
    pub(crate) const fn from_unit(v: Vec3) -> Self {
        Self { v }
    }

    /// From any non-zero finite vector (normalised).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a zero or non-finite vector.
    pub fn from_vec3(v: Vec3) -> Result<Self> {
        if !v.is_finite() {
            return Err(invalid("v", format!("must be finite, got {v:?}")));
        }
        v.try_normalize()
            .map(Self::from_unit)
            .ok_or_else(|| invalid("v", "must be non-zero"))
    }

    /// From the off-boresight angle `θ` and azimuth `φ` (radians).
    ///
    /// ```
    /// use antenna_r::Direction;
    /// let d = Direction::from_theta_phi(0.5, 1.0);
    /// assert!((d.theta_rad() - 0.5).abs() < 1e-12);
    /// assert!((d.phi_rad() - 1.0).abs() < 1e-12);
    /// ```
    #[inline]
    #[must_use]
    pub fn from_theta_phi(theta_rad: f64, phi_rad: f64) -> Self {
        let (st, ct) = theta_rad.sin_cos();
        let (sp, cp) = phi_rad.sin_cos();
        Self::from_unit(Vec3::new(st * cp, st * sp, ct))
    }

    /// From direction cosines in the front hemisphere (`w = +√(1 − u² − v²)`).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// when `u² + v² > 1` or a value is not finite.
    pub fn from_uv(u: f64, v: f64) -> Result<Self> {
        Ok(Self::from_unit(Vec3::new(u, v, w_from_uv(u, v)?)))
    }

    /// From direction cosines in the back hemisphere (`w = −√(1 − u² − v²)`).
    ///
    /// # Errors
    /// Same as [`Direction::from_uv`].
    pub fn from_uv_back(u: f64, v: f64) -> Result<Self> {
        Ok(Self::from_unit(Vec3::new(u, v, -w_from_uv(u, v)?)))
    }

    /// From X/Y tilts `atan2(u, w)` and `atan2(v, w)`, front hemisphere.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// when a tilt is not in `(-π/2, π/2)`.
    pub fn from_tilts(tilt_x_rad: f64, tilt_y_rad: f64) -> Result<Self> {
        if !(tilt_x_rad.abs() < FRAC_PI_2 && tilt_y_rad.abs() < FRAC_PI_2) {
            return Err(invalid(
                "tilt",
                format!("tilts must be in (-90°, 90°), got ({tilt_x_rad}, {tilt_y_rad}) rad"),
            ));
        }
        Self::from_vec3(Vec3::new(tilt_x_rad.tan(), tilt_y_rad.tan(), 1.0))
    }

    /// TICRA "elevation over azimuth" angles (grid type 4):
    /// `r̂ = (−sin Az cos El, sin El, cos Az cos El)`.
    #[must_use]
    pub fn from_el_over_az(az_rad: f64, el_rad: f64) -> Self {
        let (sa, ca) = az_rad.sin_cos();
        let (se, ce) = el_rad.sin_cos();
        Self::from_unit(Vec3::new(-sa * ce, se, ca * ce))
    }

    /// TICRA "elevation and azimuth" angles (grid type 5):
    /// `Az = −θ cos φ`, `El = θ sin φ`.
    #[must_use]
    pub fn from_el_and_az(az_rad: f64, el_rad: f64) -> Self {
        let theta = az_rad.hypot(el_rad);
        let phi = el_rad.atan2(-az_rad);
        Self::from_theta_phi(theta, phi)
    }

    /// TICRA "azimuth over elevation" angles (grid type 6):
    /// `r̂ = (−sin Az, cos Az sin El, cos Az cos El)`.
    #[must_use]
    pub fn from_az_over_el(az_rad: f64, el_rad: f64) -> Self {
        let (sa, ca) = az_rad.sin_cos();
        let (se, ce) = el_rad.sin_cos();
        Self::from_unit(Vec3::new(-sa, ca * se, ca * ce))
    }

    /// Unit vector in the antenna frame.
    #[inline]
    #[must_use]
    pub const fn vec3(self) -> Vec3 {
        self.v
    }

    /// Direction cosine along `+X`.
    #[inline]
    #[must_use]
    pub const fn u(self) -> f64 {
        self.v.x
    }

    /// Direction cosine along `+Y`.
    #[inline]
    #[must_use]
    pub const fn v(self) -> f64 {
        self.v.y
    }

    /// Direction cosine along `+Z` (`cos θ`).
    #[inline]
    #[must_use]
    pub const fn w(self) -> f64 {
        self.v.z
    }

    /// Off-boresight angle `θ` in `[0, π]`.
    #[inline]
    #[must_use]
    pub fn theta_rad(self) -> f64 {
        self.v.x.hypot(self.v.y).atan2(self.v.z)
    }

    /// Azimuth `φ` in `[0, 2π)` (0 on the boresight).
    #[inline]
    #[must_use]
    pub fn phi_rad(self) -> f64 {
        wrap_two_pi(self.v.y.atan2(self.v.x))
    }

    /// TICRA "elevation over azimuth" angles `(Az, El)`, inverse of
    /// [`Direction::from_el_over_az`]. `Az ∈ (−π, π]`, `El ∈ [−π/2, π/2]`.
    #[must_use]
    pub fn to_el_over_az(self) -> (f64, f64) {
        ((-self.v.x).atan2(self.v.z), self.v.y.clamp(-1.0, 1.0).asin())
    }

    /// TICRA "elevation and azimuth" angles `(Az, El)`, inverse of
    /// [`Direction::from_el_and_az`].
    #[must_use]
    pub fn to_el_and_az(self) -> (f64, f64) {
        let theta = self.theta_rad();
        let (sp, cp) = self.v.y.atan2(self.v.x).sin_cos();
        (-theta * cp, theta * sp)
    }

    /// TICRA "azimuth over elevation" angles `(Az, El)`, inverse of
    /// [`Direction::from_az_over_el`]. `Az ∈ [−π/2, π/2]`, `El ∈ (−π, π]`.
    #[must_use]
    pub fn to_az_over_el(self) -> (f64, f64) {
        ((-self.v.x).clamp(-1.0, 1.0).asin(), self.v.y.atan2(self.v.z))
    }

    /// Angle between two directions in `[0, π]` (robust `atan2` form).
    #[inline]
    #[must_use]
    pub fn angle_to(self, other: Self) -> f64 {
        self.v.angle_to(other.v)
    }

    /// The same direction in [`ScanAngles`] form.
    #[must_use]
    pub fn to_scan_angles(self) -> ScanAngles {
        ScanAngles::from_direction(self.v)
    }
}

impl From<ScanAngles> for Direction {
    fn from(s: ScanAngles) -> Self {
        Self::from_unit(s.direction())
    }
}

fn w_from_uv(u: f64, v: f64) -> Result<f64> {
    let w2 = 1.0 - u.mul_add(u, v * v);
    if !(u.is_finite() && v.is_finite()) || w2 < -UV_TOLERANCE {
        return Err(invalid(
            "uv",
            format!("need finite u, v with u²+v² ≤ 1, got ({u}, {v})"),
        ));
    }
    Ok(w2.max(0.0).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;
    use std::f64::consts::PI;

    fn assert_same(a: Direction, b: Direction) {
        assert!(a.angle_to(b) < 1e-12, "{a:?} != {b:?}");
    }

    #[test]
    fn theta_phi_round_trip_and_cosines() {
        let d = Direction::from_theta_phi(0.3, 4.0);
        assert_close(d.theta_rad(), 0.3, 1e-14);
        assert_close(d.phi_rad(), 4.0, 1e-14);
        assert_close(d.u(), 0.3_f64.sin() * 4.0_f64.cos(), 1e-15);
        assert_close(d.v(), 0.3_f64.sin() * 4.0_f64.sin(), 1e-15);
        assert_close(d.w(), 0.3_f64.cos(), 1e-15);
        assert_close(d.vec3().norm(), 1.0, 1e-15);
    }

    #[test]
    fn boresight_and_back_direction() {
        assert_eq!(Direction::BORESIGHT.theta_rad(), 0.0);
        assert_eq!(Direction::BORESIGHT.phi_rad(), 0.0);
        assert_close(Direction::from_theta_phi(PI, 0.0).theta_rad(), PI, 1e-15);
    }

    #[test]
    fn non_finite_angles_propagate_nan() {
        assert!(Direction::from_theta_phi(f64::NAN, 0.0).theta_rad().is_nan());
    }

    #[test]
    fn from_vec3_normalises_and_rejects_bad_vectors() {
        let d = Direction::from_vec3(Vec3::new(0.0, 0.0, 5.0)).unwrap();
        assert_same(d, Direction::BORESIGHT);
        assert!(matches!(
            Direction::from_vec3(Vec3::ZERO),
            Err(AntennaError::InvalidParameter { name: "v", .. })
        ));
        assert!(matches!(
            Direction::from_vec3(Vec3::new(f64::NAN, 0.0, 1.0)),
            Err(AntennaError::InvalidParameter { name: "v", .. })
        ));
    }

    #[test]
    fn uv_front_and_back_hemispheres() {
        let f = Direction::from_uv(0.6, 0.0).unwrap();
        assert_close(f.w(), 0.8, 1e-15);
        let b = Direction::from_uv_back(0.6, 0.0).unwrap();
        assert_close(b.w(), -0.8, 1e-15);
        let horizon = Direction::from_uv(1.0, 0.0).unwrap();
        assert_close(horizon.theta_rad(), FRAC_PI_2, 1e-15);
    }

    #[test]
    fn uv_outside_unit_circle_or_non_finite_is_rejected() {
        for (u, v) in [(0.8, 0.8), (f64::NAN, 0.0), (0.0, f64::INFINITY)] {
            assert!(matches!(
                Direction::from_uv(u, v),
                Err(AntennaError::InvalidParameter { name: "uv", .. })
            ));
        }
        assert!(Direction::from_uv_back(2.0, 0.0).is_err());
    }

    #[test]
    fn tilts_match_scan_angle_convention() {
        let d = Direction::from_tilts(0.2, -0.1).unwrap();
        let s = d.to_scan_angles();
        assert_close(s.tilt_x_rad, 0.2, 1e-14);
        assert_close(s.tilt_y_rad, -0.1, 1e-14);
        assert!(Direction::from_tilts(FRAC_PI_2, 0.0).is_err());
        assert!(Direction::from_tilts(0.0, f64::NAN).is_err());
    }

    #[test]
    fn el_over_az_matches_ticra_definition_and_round_trips() {
        let d = Direction::from_el_over_az(0.4, 0.2);
        assert_close(d.u(), -(0.4_f64.sin()) * 0.2_f64.cos(), 1e-15);
        assert_close(d.v(), 0.2_f64.sin(), 1e-15);
        let (az, el) = d.to_el_over_az();
        assert_close(az, 0.4, 1e-14);
        assert_close(el, 0.2, 1e-14);
    }

    #[test]
    fn az_over_el_matches_ticra_definition_and_round_trips() {
        let d = Direction::from_az_over_el(-0.3, 0.5);
        assert_close(d.u(), 0.3_f64.sin(), 1e-15);
        assert_close(d.v(), 0.3_f64.cos() * 0.5_f64.sin(), 1e-15);
        let (az, el) = d.to_az_over_el();
        assert_close(az, -0.3, 1e-14);
        assert_close(el, 0.5, 1e-14);
    }

    #[test]
    fn el_and_az_matches_ticra_definition_and_round_trips() {
        let theta = 0.25;
        let phi = 1.1;
        let d = Direction::from_theta_phi(theta, phi);
        let (az, el) = d.to_el_and_az();
        assert_close(az, -theta * phi.cos(), 1e-14);
        assert_close(el, theta * phi.sin(), 1e-14);
        assert_same(Direction::from_el_and_az(az, el), d);
    }

    #[test]
    fn scan_angles_round_trip() {
        let d = Direction::from_theta_phi(1.2, 2.5);
        let back: Direction = d.to_scan_angles().into();
        assert_same(back, d);
    }

    #[test]
    fn angle_between_directions() {
        let a = Direction::from_theta_phi(0.1, 0.0);
        let b = Direction::from_theta_phi(0.1, PI);
        assert_close(a.angle_to(b), 0.2, 1e-14);
    }
}
