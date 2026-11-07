//! Grating-lobe prediction for periodic planar lattices.
//!
//! For a lattice with primitive vectors `a1`, `a2` (X-Y plane, metres) the
//! array factor repeats at `(u, v) = (u0, v0) + λ (m b1 + n b2)`, where
//! `a_i · b_j = δ_ij` (reciprocal lattice) and `(u0, v0)` is the scan
//! direction. A repeat with `u² + v² ≤ 1` is a visible grating lobe.

use crate::direction::Direction;
use crate::error::{Result, invalid, require_positive};

/// Positions `(u, v)` of the visible grating lobes (main beam excluded).
///
/// ```
/// use antenna_r::{Direction, array::grating::grating_lobes};
/// // One-wavelength square lattice at broadside: four grating lobes at the horizon.
/// let lobes = grating_lobes((1.0, 0.0), (0.0, 1.0), 1.0, Direction::BORESIGHT).unwrap();
/// assert_eq!(lobes.len(), 4);
/// ```
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for a non-positive wavelength or degenerate (collinear or non-finite)
/// lattice vectors.
pub fn grating_lobes(a1: (f64, f64), a2: (f64, f64), wavelength_m: f64, scan: Direction) -> Result<Vec<(f64, f64)>> {
    let lambda = require_positive("wavelength_m", wavelength_m)?;
    let det = a1.0 * a2.1 - a1.1 * a2.0;
    if !(det.is_finite() && det.abs() > 0.0) {
        return Err(invalid("a1", "lattice vectors must be finite and not collinear"));
    }
    let b1 = (a2.1 / det, -a2.0 / det);
    let b2 = (-a1.1 / det, a1.0 / det);
    // |m| ≤ |Δu · a1| / λ ≤ 2 |a1| / λ, likewise for n.
    let bound = |a: (f64, f64)| (2.0 * a.0.hypot(a.1) / lambda).floor() as i64 + 1;
    let (m_max, n_max) = (bound(a1), bound(a2));
    let (u0, v0) = (scan.u(), scan.v());
    let mut lobes = Vec::new();
    for m in -m_max..=m_max {
        for n in -n_max..=n_max {
            if m == 0 && n == 0 {
                continue;
            }
            let (mf, nf) = (m as f64, n as f64);
            let u = u0 + lambda * (mf * b1.0 + nf * b2.0);
            let v = v0 + lambda * (mf * b1.1 + nf * b2.1);
            if u * u + v * v <= 1.0 + 1e-12 {
                lobes.push((u, v));
            }
        }
    }
    Ok(lobes)
}

/// Whether any grating lobe is visible (see [`grating_lobes`]).
///
/// # Errors
/// As [`grating_lobes`].
pub fn has_grating_lobe(a1: (f64, f64), a2: (f64, f64), wavelength_m: f64, scan: Direction) -> Result<bool> {
    Ok(!grating_lobes(a1, a2, wavelength_m, scan)?.is_empty())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    #[test]
    fn half_wavelength_square_lattice_has_none_up_to_horizon() {
        // d = λ/2 is grating-lobe free for any scan in the front hemisphere
        // except exactly at the horizon (u0 = 1 puts a lobe at u = −1).
        let scan = Direction::from_theta_phi(60f64.to_radians(), 0.0);
        assert!(!has_grating_lobe((0.5, 0.0), (0.0, 0.5), 1.0, scan).unwrap());
    }

    #[test]
    fn linear_spacing_rule_d_over_lambda_one_over_one_plus_sin() {
        // Grating lobe appears when d > λ / (1 + sin θ0).
        let theta = 30f64.to_radians();
        let scan = Direction::from_theta_phi(theta, 0.0);
        let limit = 1.0 / (1.0 + theta.sin());
        assert!(!has_grating_lobe((limit * 0.99, 0.0), (0.0, 0.5), 1.0, scan).unwrap());
        let lobes = grating_lobes((limit * 1.01, 0.0), (0.0, 0.5), 1.0, scan).unwrap();
        assert_eq!(lobes.len(), 1);
        assert_close(lobes[0].0, theta.sin() - 1.0 / (limit * 1.01), 1e-12);
    }

    #[test]
    fn triangular_lattice_reciprocal_vectors() {
        let d = 1.0;
        let a1 = (d, 0.0);
        let a2 = (0.5 * d, 3f64.sqrt() / 2.0 * d);
        // Equilateral lattice with spacing λ: six lobes at |Δu| = 2/√3 > 1, none visible at broadside.
        assert!(!has_grating_lobe(a1, a2, 1.0, Direction::BORESIGHT).unwrap());
        // With spacing 1.2 λ the six first-ring lobes lie at 2/(√3·1.2) ≈ 0.96 < 1.
        let a1 = (1.2, 0.0);
        let a2 = (0.6, 3f64.sqrt() / 2.0 * 1.2);
        assert_eq!(grating_lobes(a1, a2, 1.0, Direction::BORESIGHT).unwrap().len(), 6);
    }

    #[test]
    fn invalid_lattices_are_rejected() {
        assert!(matches!(
            grating_lobes((1.0, 0.0), (2.0, 0.0), 1.0, Direction::BORESIGHT),
            Err(AntennaError::InvalidParameter { name: "a1", .. })
        ));
        assert!(grating_lobes((f64::NAN, 0.0), (0.0, 1.0), 1.0, Direction::BORESIGHT).is_err());
        assert!(grating_lobes((1.0, 0.0), (0.0, 1.0), 0.0, Direction::BORESIGHT).is_err());
        assert!(has_grating_lobe((1.0, 0.0), (0.0, 1.0), -1.0, Direction::BORESIGHT).is_err());
    }
}
