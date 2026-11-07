//! Gain towards ECEF targets using the frames of `geometry_r`.
//!
//! - Satellite antennas: the platform body frame comes from a
//!   [`SatelliteFrame`]; the antenna is mounted with a fixed rotation whose
//!   columns are the antenna axes expressed in the body frame
//!   ([`Mat3::IDENTITY`] when the antenna frame is the body frame).
//! - Ground antennas: the mount frame is the local ENU frame of a
//!   [`LocalFrame`], matching [`crate::steer::Mount::AzEl`].

use geometry_r::coords::{Ecef, LocalFrame};
use geometry_r::math::Mat3;
use geometry_r::satellite::SatelliteFrame;

use crate::direction::Direction;
use crate::error::{Result, invalid};
use crate::pattern::GainPattern;

/// Antenna-frame direction from a satellite to an ECEF target.
///
/// # Errors
/// [`AntennaError::Geometry`](crate::AntennaError::Geometry) when the target
/// coincides with the satellite.
pub fn direction_from_satellite(frame: &SatelliteFrame, mounting: &Mat3, target: Ecef) -> Result<Direction> {
    let body = frame.direction_to(target)?;
    Ok(Direction::from_unit(mounting.mul_vec_t(body)))
}

/// Gain (dBi) of a satellite antenna towards an ECEF target.
///
/// ```
/// use antenna_r::{GainPattern, frames::gain_from_satellite, parabolic::GaussianBeam};
/// use geometry_r::{math::Mat3, satellite::{SatelliteFrame, NadirReference, geo_satellite_position, sub_satellite_point}};
/// let sat = geo_satellite_position(0.0);
/// let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric).unwrap();
/// let beam = GaussianBeam::new(40.0, 1f64.to_radians()).unwrap();
/// let nadir = sub_satellite_point(sat).to_ecef();
/// let g = gain_from_satellite(&beam, &frame, &Mat3::IDENTITY, nadir).unwrap();
/// assert!((g - 40.0).abs() < 1e-6);
/// ```
///
/// # Errors
/// As [`direction_from_satellite`].
pub fn gain_from_satellite<P: GainPattern + ?Sized>(
    pattern: &P,
    frame: &SatelliteFrame,
    mounting: &Mat3,
    target: Ecef,
) -> Result<f64> {
    Ok(pattern.gain_dbi(direction_from_satellite(frame, mounting, target)?))
}

/// Gains towards many targets; a target at the satellite position gives `NaN`.
///
/// # Panics
/// Panics if `targets.len() != out.len()`.
pub fn gains_from_satellite<P: GainPattern + ?Sized>(
    pattern: &P,
    frame: &SatelliteFrame,
    mounting: &Mat3,
    targets: &[Ecef],
    out: &mut [f64],
) {
    assert_eq!(targets.len(), out.len(), "targets and out must have the same length");
    for (o, t) in out.iter_mut().zip(targets) {
        *o = gain_from_satellite(pattern, frame, mounting, *t).unwrap_or(f64::NAN);
    }
}

/// ENU (mount-frame) direction from a ground station to an ECEF target.
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// when the target coincides with the station.
pub fn direction_from_station(station: &LocalFrame, target: Ecef) -> Result<Direction> {
    Direction::from_vec3(station.to_enu(target)).map_err(|_| invalid("target", "coincides with the station"))
}

/// Gain (dBi) of a ground antenna, whose pattern is expressed in the ENU
/// frame (typically a [`Steered`](crate::steer::Steered) with an az/el
/// mount), towards an ECEF target.
///
/// # Errors
/// As [`direction_from_station`].
pub fn gain_from_station<P: GainPattern + ?Sized>(pattern: &P, station: &LocalFrame, target: Ecef) -> Result<f64> {
    Ok(pattern.gain_dbi(direction_from_station(station, target)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;
    use crate::parabolic::GaussianBeam;
    use crate::steer::{Mount, Steered};
    use geometry_r::coords::Geodetic;
    use geometry_r::math::Vec3;
    use geometry_r::satellite::{NadirReference, geo_satellite_position, sub_satellite_point};

    fn beam() -> GaussianBeam {
        GaussianBeam::new(40.0, 1f64.to_radians()).unwrap()
    }

    #[test]
    fn satellite_beam_peaks_at_nadir_and_mounting_rotates_it() {
        let sat = geo_satellite_position(0.0);
        let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric).unwrap();
        let nadir = sub_satellite_point(sat).to_ecef();
        assert_close(
            gain_from_satellite(&beam(), &frame, &Mat3::IDENTITY, nadir).unwrap(),
            40.0,
            1e-6,
        );
        // Antenna boresight tilted 0.5° about body X: nadir is now 0.5° off.
        let tilt = Mat3::rot_x(0.5f64.to_radians());
        let g = gain_from_satellite(&beam(), &frame, &tilt, nadir).unwrap();
        assert_close(g, 37.0, 1e-6);
        let d = direction_from_satellite(&frame, &tilt, nadir).unwrap();
        assert_close(d.theta_rad(), 0.5f64.to_radians(), 1e-9);
    }

    #[test]
    fn satellite_errors_and_batch_nan() {
        let sat = geo_satellite_position(0.0);
        let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric).unwrap();
        assert!(matches!(
            gain_from_satellite(&beam(), &frame, &Mat3::IDENTITY, sat),
            Err(AntennaError::Geometry(_))
        ));
        let nadir = sub_satellite_point(sat).to_ecef();
        let mut out = [0.0; 2];
        gains_from_satellite(&beam(), &frame, &Mat3::IDENTITY, &[nadir, sat], &mut out);
        assert_close(out[0], 40.0, 1e-6);
        assert!(out[1].is_nan());
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn satellite_batch_panics_on_length_mismatch() {
        let sat = geo_satellite_position(0.0);
        let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric).unwrap();
        gains_from_satellite(&beam(), &frame, &Mat3::IDENTITY, &[sat], &mut []);
    }

    #[test]
    fn station_antenna_tracks_a_satellite_and_sees_a_neighbour_off_axis() {
        let station = LocalFrame::new(&Geodetic::from_degrees(45.0, 2.0, 0.0).unwrap());
        let wanted = geo_satellite_position(5f64.to_radians());
        let neighbour = geo_satellite_position(7f64.to_radians());
        let mut dish = Steered::new(beam(), Mount::AzEl);
        dish.point_at(direction_from_station(&station, wanted).unwrap())
            .unwrap();
        assert_close(gain_from_station(&dish, &station, wanted).unwrap(), 40.0, 1e-9);
        let off = gain_from_station(&dish, &station, neighbour).unwrap();
        let sep = direction_from_station(&station, wanted)
            .unwrap()
            .angle_to(direction_from_station(&station, neighbour).unwrap());
        assert_close(off, 40.0 - 12.0 * (sep / 1f64.to_radians()).powi(2), 1e-9);
        assert!(matches!(
            gain_from_station(&dish, &station, station.origin),
            Err(AntennaError::InvalidParameter { name: "target", .. })
        ));
        let up = station.from_enu(Vec3::new(0.0, 0.0, 1000.0));
        assert_close(direction_from_station(&station, up).unwrap().w(), 1.0, 1e-9);
    }
}
