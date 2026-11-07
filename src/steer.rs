//! Steerable antennas: any [`GainPattern`] pointed by a mount, with axis
//! limits, roll, pointing error and an optional electronic scan loss.
//!
//! Two frames are involved. The **mount frame** is fixed to the platform
//! (for a ground az/el mount: X east, Y north, Z up, i.e. ENU). The
//! **antenna frame** moves with the antenna (boresight `+Z`). Gains are
//! requested in the mount frame and evaluated in the antenna frame:
//! `d_ant = Rᵀ d_mount`, where the columns of `R` are the antenna axes
//! expressed in the mount frame.

use std::f64::consts::{FRAC_PI_2, PI};

use geometry_r::math::{Mat3, Vec3, wrap_two_pi};

use crate::direction::Direction;
use crate::error::{AntennaError, Result, invalid, require_finite};
use crate::pattern::GainPattern;

/// How the antenna is moved. Each mount has two axis angles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mount {
    /// Free rotation (hexapod, electronic steering): axes `(θ, φ)` of the
    /// boresight in the mount frame, reached by the smallest rotation from
    /// `+Z`.
    Free,
    /// Azimuth/elevation: azimuth clockwise from `+Y` towards `+X` in
    /// `[0, 2π)`, elevation from the X-Y plane towards `+Z`. With an ENU
    /// mount frame these are the usual azimuth from north and elevation.
    /// Antenna `+X` points towards increasing elevation. At the zenith
    /// (keyhole) the azimuth is undefined and taken as 0.
    AzEl,
    /// X-over-Y: rotation `x` about the mount X axis, then `y` about the
    /// rotated Y axis; boresight `(sin y, −sin x cos y, cos x cos y)`.
    XOverY,
}

impl Mount {
    /// Rotation whose columns are the antenna axes in the mount frame for
    /// axis angles `(a1, a2)`.
    #[must_use]
    pub fn rotation(self, a1: f64, a2: f64) -> Mat3 {
        match self {
            Self::Free => Mat3::rot_z(a2).mul_mat(&Mat3::rot_y(a1)).mul_mat(&Mat3::rot_z(-a2)),
            Self::AzEl => {
                let (sa, ca) = a1.sin_cos();
                let (se, ce) = a2.sin_cos();
                let z = Vec3::new(ce * sa, ce * ca, se);
                let x = Vec3::new(-se * sa, -se * ca, ce);
                Mat3::from_cols(x, z.cross(x), z)
            }
            Self::XOverY => Mat3::rot_x(a1).mul_mat(&Mat3::rot_y(a2)),
        }
    }

    /// Axis angles that point the boresight at `target` (mount frame).
    #[must_use]
    pub fn axes_for(self, target: Direction) -> (f64, f64) {
        let b = target.vec3();
        match self {
            Self::Free => (target.theta_rad(), target.phi_rad()),
            Self::AzEl => (wrap_two_pi(b.x.atan2(b.y)), b.z.clamp(-1.0, 1.0).asin()),
            Self::XOverY => ((-b.y).atan2(b.z), b.x.clamp(-1.0, 1.0).asin()),
        }
    }
}

/// Allowed range of one axis (radians, inclusive).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisRange {
    /// Lowest angle.
    pub min_rad: f64,
    /// Highest angle.
    pub max_rad: f64,
}

impl AxisRange {
    /// Validated range.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for non-finite or reversed limits.
    pub fn new(min_rad: f64, max_rad: f64) -> Result<Self> {
        require_finite("min_rad", min_rad)?;
        require_finite("max_rad", max_rad)?;
        if min_rad > max_rad {
            return Err(invalid("max_rad", format!("{max_rad} < min_rad {min_rad}")));
        }
        Ok(Self { min_rad, max_rad })
    }

    fn contains(&self, a: f64) -> bool {
        a >= self.min_rad && a <= self.max_rad
    }

    fn clamp(&self, a: f64) -> f64 {
        a.clamp(self.min_rad, self.max_rad)
    }
}

/// Steering limits: optional per-axis ranges and a maximum scan angle (the
/// angle between the boresight and the mount `+Z`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteeringLimits {
    /// Range of the first axis.
    pub axis1: Option<AxisRange>,
    /// Range of the second axis.
    pub axis2: Option<AxisRange>,
    /// Largest scan angle in `[0, π]`.
    pub max_scan_rad: f64,
}

impl Default for SteeringLimits {
    fn default() -> Self {
        Self {
            axis1: None,
            axis2: None,
            max_scan_rad: PI,
        }
    }
}

/// Electronic scan loss: gain changes by `10 n log10(cos θs)` and, if
/// `broaden` is set, the beam widens by `1 / cos θs` in the scan plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanLoss {
    exponent: f64,
    broaden: bool,
}

impl ScanLoss {
    /// Scan loss with exponent `n ≥ 0` (1 for projected-aperture loss).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for a negative or non-finite exponent.
    pub fn new(exponent: f64, broaden: bool) -> Result<Self> {
        if !(exponent.is_finite() && exponent >= 0.0) {
            return Err(invalid("exponent", format!("must be finite and ≥ 0, got {exponent}")));
        }
        Ok(Self { exponent, broaden })
    }

    /// Loss term in dB (≤ 0) for a scan angle below 90°.
    #[must_use]
    pub fn gain_change_db(&self, scan_rad: f64) -> f64 {
        self.exponent * rf_r::db::to_db(scan_rad.cos())
    }
}

/// Precomputed scan-plane compression for beam broadening.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Broadening {
    e_scan: Vec3,
    e_cross: Vec3,
    cos_scan: f64,
}

/// A pattern carried by a mount.
#[derive(Debug, Clone, PartialEq)]
pub struct Steered<P> {
    pattern: P,
    mount: Mount,
    limits: SteeringLimits,
    scan_loss: Option<ScanLoss>,
    axes: (f64, f64),
    roll_rad: f64,
    bias: Mat3,
    rot: Mat3,
    scan_rad: f64,
    scan_change_db: f64,
    broadening: Option<Broadening>,
}

impl<P: GainPattern> Steered<P> {
    /// Mounts `pattern`, initially pointed along the mount `+Z`, with no
    /// limits, roll, pointing error or scan loss.
    #[must_use]
    pub fn new(pattern: P, mount: Mount) -> Self {
        let mut s = Self {
            pattern,
            mount,
            limits: SteeringLimits::default(),
            scan_loss: None,
            axes: mount.axes_for(Direction::BORESIGHT),
            roll_rad: 0.0,
            bias: Mat3::IDENTITY,
            rot: Mat3::IDENTITY,
            scan_rad: 0.0,
            scan_change_db: 0.0,
            broadening: None,
        };
        s.update();
        s
    }

    /// Sets the steering limits (checked by later pointing calls).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for a maximum scan outside `[0, π]`.
    pub fn with_limits(mut self, limits: SteeringLimits) -> Result<Self> {
        if !(0.0..=PI).contains(&limits.max_scan_rad) {
            return Err(invalid(
                "max_scan_rad",
                format!("must be in [0, π], got {}", limits.max_scan_rad),
            ));
        }
        self.limits = limits;
        Ok(self)
    }

    /// Enables an electronic scan loss. The maximum scan is reduced below
    /// 90° if needed (the loss is infinite at 90°).
    #[must_use]
    pub fn with_scan_loss(mut self, scan_loss: ScanLoss) -> Self {
        self.scan_loss = Some(scan_loss);
        self.limits.max_scan_rad = self.limits.max_scan_rad.min(FRAC_PI_2 - 1e-9);
        self.update();
        self
    }

    /// Points the boresight at `target` (mount frame).
    ///
    /// # Errors
    /// [`AntennaError::OutOfRange`] when the required axes or scan angle
    /// violate the limits; the pointing is then unchanged.
    pub fn point_at(&mut self, target: Direction) -> Result<()> {
        let (a1, a2) = self.mount.axes_for(target);
        self.set_axes(a1, a2)
    }

    /// Sets the axis angles directly.
    ///
    /// # Errors
    /// [`AntennaError::OutOfRange`] when an axis or the scan angle violates
    /// the limits, [`AntennaError::InvalidParameter`] for non-finite angles.
    pub fn set_axes(&mut self, a1: f64, a2: f64) -> Result<()> {
        require_finite("a1", a1)?;
        require_finite("a2", a2)?;
        let check = |r: Option<AxisRange>, a: f64, name: &str| match r {
            Some(r) if !r.contains(a) => Err(AntennaError::OutOfRange(format!(
                "{name} = {a} rad outside [{}, {}]",
                r.min_rad, r.max_rad
            ))),
            _ => Ok(()),
        };
        check(self.limits.axis1, a1, "axis 1")?;
        check(self.limits.axis2, a2, "axis 2")?;
        let scan = self.mount.rotation(a1, a2).mul_vec(Vec3::Z).angle_to(Vec3::Z);
        if scan > self.limits.max_scan_rad {
            return Err(AntennaError::OutOfRange(format!(
                "scan angle {scan} rad above {} rad",
                self.limits.max_scan_rad
            )));
        }
        self.axes = (a1, a2);
        self.update();
        Ok(())
    }

    /// Points as close to `target` as the limits allow: the scan angle is
    /// pulled back to the cone, then each axis is clamped to its range.
    pub fn point_at_clamped(&mut self, target: Direction) {
        let max = self.limits.max_scan_rad;
        let target = if target.theta_rad() > max {
            Direction::from_theta_phi(max, target.phi_rad())
        } else {
            target
        };
        let (a1, a2) = self.mount.axes_for(target);
        let clamp = |r: Option<AxisRange>, a: f64| r.map_or(a, |r| r.clamp(a));
        self.axes = (clamp(self.limits.axis1, a1), clamp(self.limits.axis2, a2));
        self.update();
    }

    /// Rotation of the antenna about its own boresight (e.g. polarisation).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for a non-finite angle.
    pub fn set_roll(&mut self, roll_rad: f64) -> Result<()> {
        self.roll_rad = require_finite("roll_rad", roll_rad)?;
        self.update();
        Ok(())
    }

    /// Fixed pointing error: the boresight is tilted by `error_x_rad`
    /// towards antenna `+X` and `error_y_rad` towards antenna `+Y`.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`] for non-finite angles.
    pub fn set_pointing_error(&mut self, error_x_rad: f64, error_y_rad: f64) -> Result<()> {
        require_finite("error_x_rad", error_x_rad)?;
        require_finite("error_y_rad", error_y_rad)?;
        self.bias = Mat3::rot_y(error_x_rad).mul_mat(&Mat3::rot_x(-error_y_rad));
        self.update();
        Ok(())
    }

    /// The carried pattern.
    #[must_use]
    pub fn pattern(&self) -> &P {
        &self.pattern
    }

    /// Current axis angles.
    #[must_use]
    pub fn axes(&self) -> (f64, f64) {
        self.axes
    }

    /// Actual boresight in the mount frame (pointing error included).
    #[must_use]
    pub fn boresight(&self) -> Direction {
        Direction::from_unit(self.rot.mul_vec(Vec3::Z))
    }

    /// Angle between the commanded boresight and the mount `+Z`.
    #[must_use]
    pub fn scan_angle_rad(&self) -> f64 {
        self.scan_rad
    }

    /// Antenna-frame direction of a mount-frame direction.
    #[inline]
    #[must_use]
    pub fn to_antenna_frame(&self, dir: Direction) -> Direction {
        Direction::from_unit(self.rot.mul_vec_t(dir.vec3()))
    }

    fn update(&mut self) {
        let mount = self.mount.rotation(self.axes.0, self.axes.1);
        self.scan_rad = mount.mul_vec(Vec3::Z).angle_to(Vec3::Z);
        self.rot = mount.mul_mat(&self.bias).mul_mat(&Mat3::rot_z(self.roll_rad));
        self.scan_change_db = self.scan_loss.map_or(0.0, |s| s.gain_change_db(self.scan_rad));
        self.broadening = self.scan_loss.filter(|s| s.broaden).and_then(|_| {
            // Mount +Z seen from the antenna: its transverse part is the scan plane.
            let z = self.rot.mul_vec_t(Vec3::Z);
            let e_scan = Vec3::new(z.x, z.y, 0.0).try_normalize()?;
            Some(Broadening {
                e_scan,
                e_cross: Vec3::Z.cross(e_scan),
                cos_scan: self.scan_rad.cos(),
            })
        });
    }

    fn broaden(&self, d: Vec3) -> Vec3 {
        let Some(b) = self.broadening else { return d };
        let a = d.dot(b.e_scan) * b.cos_scan;
        let c = d.dot(b.e_cross);
        let w = (1.0 - a * a - c * c).max(0.0).sqrt().copysign(d.z);
        b.e_scan * a + b.e_cross * c + Vec3::Z * w
    }
}

impl<P: GainPattern> GainPattern for Steered<P> {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let d = self.broaden(self.rot.mul_vec_t(dir.vec3()));
        self.pattern.gain_dbi(Direction::from_unit(d)) + self.scan_change_db
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.pattern.peak_gain_dbi() + self.scan_change_db
    }
}

/// Gain lost by mispointing a pattern by `error_rad` in the plane `phi_rad`:
/// boresight gain minus the gain at `(θ = error, φ)`, in dB (≥ 0 for a beam
/// peaked on boresight).
#[must_use]
pub fn pointing_loss_db<P: GainPattern + ?Sized>(pattern: &P, error_rad: f64, phi_rad: f64) -> f64 {
    pattern.gain_dbi(Direction::BORESIGHT) - pattern.gain_dbi(Direction::from_theta_phi(error_rad, phi_rad))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::assert_close;
    use crate::parabolic::GaussianBeam;

    fn beam() -> GaussianBeam {
        GaussianBeam::new(30.0, 0.05).unwrap()
    }

    fn assert_same(a: Direction, b: Direction) {
        assert!(a.angle_to(b) < 1e-12, "{a:?} != {b:?}");
    }

    #[test]
    fn every_mount_points_its_boresight_at_the_target() {
        let target = Direction::from_theta_phi(0.7, 2.2);
        for mount in [Mount::Free, Mount::AzEl, Mount::XOverY] {
            let (a1, a2) = mount.axes_for(target);
            let r = mount.rotation(a1, a2);
            assert_same(Direction::from_unit(r.mul_vec(Vec3::Z)), target);
            // Orthonormal, right-handed axes.
            let (x, y, z) = (r.mul_vec(Vec3::X), r.mul_vec(Vec3::Y), r.mul_vec(Vec3::Z));
            assert_close(x.cross(y).dot(z), 1.0, 1e-12);
            let mut s = Steered::new(beam(), mount);
            s.point_at(target).unwrap();
            assert_close(s.gain_dbi(target), 30.0, 1e-9);
            assert_same(s.boresight(), target);
        }
    }

    #[test]
    fn az_el_axes_use_compass_azimuth() {
        // North-east at 30° elevation in ENU.
        let el = 30f64.to_radians();
        let target =
            Direction::from_vec3(Vec3::new(el.cos() * 0.5f64.sqrt(), el.cos() * 0.5f64.sqrt(), el.sin())).unwrap();
        let (az, e) = Mount::AzEl.axes_for(target);
        assert_close(az, 45f64.to_radians(), 1e-12);
        assert_close(e, el, 1e-12);
        // Keyhole: zenith gives azimuth 0.
        assert_eq!(Mount::AzEl.axes_for(Direction::BORESIGHT), (0.0, FRAC_PI_2));
    }

    #[test]
    fn x_over_y_axes_follow_definition() {
        let (x, y) = (0.3, -0.2);
        let b = Mount::XOverY.rotation(x, y).mul_vec(Vec3::Z);
        assert_close(b.x, y.sin(), 1e-15);
        assert_close(b.y, -x.sin() * y.cos(), 1e-15);
        let d = Direction::from_unit(b);
        let (ax, ay) = Mount::XOverY.axes_for(d);
        assert_close(ax, x, 1e-12);
        assert_close(ay, y, 1e-12);
    }

    #[test]
    fn off_axis_gain_follows_the_pattern() {
        let mut s = Steered::new(beam(), Mount::Free);
        let target = Direction::from_theta_phi(0.4, 1.0);
        s.point_at(target).unwrap();
        let off = Direction::from_theta_phi(0.425, 1.0);
        assert_close(s.gain_dbi(off), 27.0, 1e-6);
        assert_same(s.to_antenna_frame(target), Direction::BORESIGHT);
        assert_close(s.scan_angle_rad(), 0.4, 1e-12);
        assert_eq!(s.peak_gain_dbi(), 30.0);
        assert_eq!(s.pattern().peak_gain_dbi(), 30.0);
    }

    #[test]
    fn roll_rotates_an_elliptical_beam() {
        let p = GaussianBeam::elliptical(30.0, 0.2, 0.05, 0.0).unwrap();
        let mut s = Steered::new(p, Mount::Free);
        let along_x = Direction::from_theta_phi(0.1, 0.0);
        assert_close(s.gain_dbi(along_x), 27.0, 1e-9);
        s.set_roll(FRAC_PI_2).unwrap();
        // Major axis now along mount +Y.
        assert_close(s.gain_dbi(Direction::from_theta_phi(0.1, FRAC_PI_2)), 27.0, 1e-9);
        assert!(s.set_roll(f64::NAN).is_err());
    }

    #[test]
    fn pointing_error_moves_the_beam() {
        let mut s = Steered::new(beam(), Mount::Free);
        s.set_pointing_error(0.025, 0.0).unwrap();
        assert_close(s.gain_dbi(Direction::BORESIGHT), 27.0, 1e-9);
        assert_same(s.boresight(), Direction::from_theta_phi(0.025, 0.0));
        s.set_pointing_error(0.0, 0.025).unwrap();
        assert_same(s.boresight(), Direction::from_theta_phi(0.025, FRAC_PI_2));
        assert!(s.set_pointing_error(f64::NAN, 0.0).is_err());
        assert!(s.set_pointing_error(0.0, f64::INFINITY).is_err());
        assert_close(pointing_loss_db(&beam(), 0.025, 0.0), 3.0, 1e-9);
    }

    #[test]
    fn limits_reject_or_clamp_pointing() {
        let limits = SteeringLimits {
            axis1: Some(AxisRange::new(0.0, PI).unwrap()),
            axis2: Some(AxisRange::new(0.1, 1.2).unwrap()),
            max_scan_rad: PI,
        };
        let mut s = Steered::new(beam(), Mount::AzEl).with_limits(limits).unwrap();
        let low = Direction::from_vec3(Vec3::new(0.0, 1.0, 0.0)).unwrap(); // el = 0
        assert!(matches!(s.point_at(low), Err(AntennaError::OutOfRange(_))));
        let west = Direction::from_vec3(Vec3::new(-1.0, 0.0, 1.0)).unwrap(); // az = 270°
        assert!(matches!(s.point_at(west), Err(AntennaError::OutOfRange(_))));
        s.point_at_clamped(low);
        assert_close(s.axes().1, 0.1, 1e-15);
        assert!(s.set_axes(f64::NAN, 0.5).is_err());
        assert!(s.set_axes(0.5, f64::NAN).is_err());
        s.set_axes(0.5, 0.5).unwrap();
        assert_eq!(s.axes(), (0.5, 0.5));
    }

    #[test]
    fn scan_cone_is_enforced() {
        let limits = SteeringLimits {
            max_scan_rad: 0.5,
            ..SteeringLimits::default()
        };
        let mut s = Steered::new(beam(), Mount::Free).with_limits(limits).unwrap();
        let far = Direction::from_theta_phi(0.8, 1.0);
        assert!(matches!(s.point_at(far), Err(AntennaError::OutOfRange(_))));
        s.point_at_clamped(far);
        assert_close(s.scan_angle_rad(), 0.5, 1e-12);
        assert_close(s.boresight().phi_rad(), 1.0, 1e-12);
        s.point_at_clamped(Direction::from_theta_phi(0.2, 1.0));
        assert_close(s.scan_angle_rad(), 0.2, 1e-12);
        let bad = SteeringLimits {
            max_scan_rad: 4.0,
            ..SteeringLimits::default()
        };
        assert!(Steered::new(beam(), Mount::Free).with_limits(bad).is_err());
    }

    #[test]
    fn axis_range_validation() {
        assert!(AxisRange::new(1.0, 0.0).is_err());
        assert!(AxisRange::new(f64::NAN, 0.0).is_err());
        assert!(AxisRange::new(0.0, f64::INFINITY).is_err());
    }

    #[test]
    fn scan_loss_reduces_gain_and_broadens_the_beam() {
        let loss = ScanLoss::new(1.0, true).unwrap();
        let mut s = Steered::new(beam(), Mount::Free).with_scan_loss(loss);
        let scan = 60f64.to_radians();
        s.point_at(Direction::from_theta_phi(scan, 0.0)).unwrap();
        let peak = 30.0 + 10.0 * 0.5f64.log10();
        assert_close(s.peak_gain_dbi(), peak, 1e-9);
        assert_close(s.gain_dbi(Direction::from_theta_phi(scan, 0.0)), peak, 1e-9);
        // In the scan plane the −3 dB point moves out to 0.025 / cos 60°.
        let edge = s.boresight().angle_to(Direction::from_theta_phi(scan + 0.05, 0.0));
        assert_close(edge, 0.05, 1e-12);
        let compressed = (0.5 * 0.05f64.sin()).asin();
        let expected = peak - 12.0 * (compressed / 0.05).powi(2);
        assert_close(s.gain_dbi(Direction::from_theta_phi(scan + 0.05, 0.0)), expected, 1e-9);
        // Across the scan plane the beam is unchanged.
        let in_plane = Vec3::new(scan.cos(), 0.0, -scan.sin());
        let cross = Direction::from_vec3(s.boresight().vec3().rotate_about(in_plane, 0.025)).unwrap();
        assert_close(s.gain_dbi(cross), peak - 3.0, 1e-9);
        // At broadside nothing changes.
        s.point_at(Direction::BORESIGHT).unwrap();
        assert_close(s.gain_dbi(Direction::from_theta_phi(0.025, 0.0)), 27.0, 1e-9);
    }

    #[test]
    fn scan_loss_without_broadening_and_limits() {
        let loss = ScanLoss::new(1.5, false).unwrap();
        assert_close(loss.gain_change_db(60f64.to_radians()), 15.0 * 0.5f64.log10(), 1e-12);
        let mut s = Steered::new(beam(), Mount::Free).with_scan_loss(loss);
        assert!(s.point_at(Direction::from_theta_phi(FRAC_PI_2, 0.0)).is_err());
        assert!(ScanLoss::new(-1.0, true).is_err());
        assert!(ScanLoss::new(f64::NAN, true).is_err());
    }
}
