//! Element lattices (positions in metres, antenna frame, centred on the
//! origin in the X-Y plane) and aperture masks.
//!
//! To give spacings in wavelengths multiply by
//! `rf_r::fspl::wavelength_m(design_freq_hz)`.

use std::f64::consts::TAU;

use geometry_r::math::Vec3;

use crate::error::{Result, invalid, require_positive};

fn require_count(name: &'static str, n: usize) -> Result<usize> {
    if n == 0 {
        Err(invalid(name, "must be at least 1"))
    } else {
        Ok(n)
    }
}

/// Moves the centroid of the points to the origin.
fn centred(mut points: Vec<Vec3>) -> Vec<Vec3> {
    let n = points.len() as f64;
    let c = points.iter().fold(Vec3::ZERO, |acc, p| acc + *p) / n;
    for p in &mut points {
        *p -= c;
    }
    points
}

/// `n` elements along X with spacing `spacing_m`.
///
/// ```
/// use antenna_r::array::lattice::linear;
/// let p = linear(3, 0.5).unwrap();
/// assert_eq!(p[0].x, -0.5);
/// ```
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for `n = 0` or a non-positive spacing.
pub fn linear(n: usize, spacing_m: f64) -> Result<Vec<Vec3>> {
    rectangular(n, 1, spacing_m, 1.0)
}

/// `nx × ny` rectangular grid, X varying fastest.
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for a zero count or a non-positive spacing.
pub fn rectangular(nx: usize, ny: usize, dx_m: f64, dy_m: f64) -> Result<Vec<Vec3>> {
    lattice(nx, ny, dx_m, dy_m, 0.0)
}

/// Triangular grid: rectangular rows with every odd row shifted by `dx/2`
/// (equilateral for `dy = dx·√3/2`).
///
/// # Errors
/// As [`rectangular`].
pub fn triangular(nx: usize, ny: usize, dx_m: f64, dy_m: f64) -> Result<Vec<Vec3>> {
    lattice(nx, ny, dx_m, dy_m, 0.5)
}

fn lattice(nx: usize, ny: usize, dx_m: f64, dy_m: f64, odd_row_shift: f64) -> Result<Vec<Vec3>> {
    require_count("nx", nx)?;
    require_count("ny", ny)?;
    require_positive("dx_m", dx_m)?;
    require_positive("dy_m", dy_m)?;
    let mut points = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        let shift = if j % 2 == 1 { odd_row_shift } else { 0.0 };
        for i in 0..nx {
            points.push(Vec3::new((i as f64 + shift) * dx_m, j as f64 * dy_m, 0.0));
        }
    }
    Ok(centred(points))
}

/// One centre element plus `n_rings` concentric rings of radius
/// `k · ring_spacing_m`, each holding `round(2π r / element_spacing_m)`
/// equally spaced elements (at least one).
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for a non-positive spacing.
pub fn concentric_rings(n_rings: usize, ring_spacing_m: f64, element_spacing_m: f64) -> Result<Vec<Vec3>> {
    require_positive("ring_spacing_m", ring_spacing_m)?;
    require_positive("element_spacing_m", element_spacing_m)?;
    let mut points = vec![Vec3::ZERO];
    for k in 1..=n_rings {
        let r = k as f64 * ring_spacing_m;
        let count = ((TAU * r / element_spacing_m).round() as usize).max(1);
        for i in 0..count {
            let (s, c) = (TAU * i as f64 / count as f64).sin_cos();
            points.push(Vec3::new(r * c, r * s, 0.0));
        }
    }
    Ok(points)
}

/// Keeps the points whose X-Y distance from the origin is at most `radius_m`.
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for a non-positive radius.
pub fn circular_mask(points: &[Vec3], radius_m: f64) -> Result<Vec<Vec3>> {
    require_positive("radius_m", radius_m)?;
    Ok(points
        .iter()
        .copied()
        .filter(|p| p.x.hypot(p.y) <= radius_m * (1.0 + 1e-12))
        .collect())
}

/// Keeps the points inside a polygon of X-Y vertices (even-odd rule).
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for fewer than 3 vertices or non-finite vertices.
pub fn polygon_mask(points: &[Vec3], polygon: &[(f64, f64)]) -> Result<Vec<Vec3>> {
    if polygon.len() < 3 {
        return Err(invalid("polygon", format!("needs ≥ 3 vertices, got {}", polygon.len())));
    }
    if polygon.iter().any(|(x, y)| !(x.is_finite() && y.is_finite())) {
        return Err(invalid("polygon", "vertices must be finite"));
    }
    Ok(points.iter().copied().filter(|p| inside(polygon, p.x, p.y)).collect())
}

fn inside(polygon: &[(f64, f64)], x: f64, y: f64) -> bool {
    let mut inside = false;
    let mut j = polygon.len() - 1;
    for (i, &(xi, yi)) in polygon.iter().enumerate() {
        let (xj, yj) = polygon[j];
        if (yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi {
            inside = !inside;
        }
        j = i;
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    #[test]
    fn linear_array_is_centred_along_x() {
        let p = linear(4, 0.5).unwrap();
        assert_eq!(p.len(), 4);
        assert_close(p[0].x, -0.75, 1e-15);
        assert_close(p[3].x, 0.75, 1e-15);
        assert!(p.iter().all(|q| q.y == 0.0 && q.z == 0.0));
    }

    #[test]
    fn rectangular_grid_is_x_fastest_and_centred() {
        let p = rectangular(3, 2, 1.0, 2.0).unwrap();
        assert_eq!(p.len(), 6);
        assert_close(p[1].x - p[0].x, 1.0, 1e-15);
        assert_close(p[3].y - p[0].y, 2.0, 1e-15);
        let c = p.iter().fold(Vec3::ZERO, |a, q| a + *q);
        assert_close(c.norm(), 0.0, 1e-12);
    }

    #[test]
    fn triangular_grid_shifts_odd_rows() {
        let p = triangular(2, 2, 1.0, 0.866).unwrap();
        assert_close(p[2].x - p[0].x, 0.5, 1e-15);
    }

    #[test]
    fn lattices_reject_invalid_sizes() {
        assert!(matches!(
            linear(0, 1.0),
            Err(AntennaError::InvalidParameter { name: "nx", .. })
        ));
        assert!(matches!(
            rectangular(1, 0, 1.0, 1.0),
            Err(AntennaError::InvalidParameter { name: "ny", .. })
        ));
        assert!(matches!(
            linear(2, 0.0),
            Err(AntennaError::InvalidParameter { name: "dx_m", .. })
        ));
        assert!(matches!(
            rectangular(2, 2, 1.0, f64::NAN),
            Err(AntennaError::InvalidParameter { name: "dy_m", .. })
        ));
    }

    #[test]
    fn rings_place_elements_on_circles() {
        let p = concentric_rings(2, 0.5, 0.5).unwrap();
        // 1 + round(2π·0.5/0.5) + round(2π·1/0.5) = 1 + 6 + 13
        assert_eq!(p.len(), 20);
        assert_close(p[1].x.hypot(p[1].y), 0.5, 1e-15);
        assert_close(p[19].x.hypot(p[19].y), 1.0, 1e-15);
        let tiny = concentric_rings(1, 0.01, 1.0).unwrap();
        assert_eq!(tiny.len(), 2);
        assert_eq!(concentric_rings(0, 1.0, 1.0).unwrap(), vec![Vec3::ZERO]);
        assert!(concentric_rings(1, 0.0, 1.0).is_err());
        assert!(concentric_rings(1, 1.0, -1.0).is_err());
    }

    #[test]
    fn circular_mask_keeps_inner_points() {
        let p = rectangular(5, 5, 1.0, 1.0).unwrap();
        let kept = circular_mask(&p, 1.0).unwrap();
        assert_eq!(kept.len(), 5);
        assert!(circular_mask(&p, 0.0).is_err());
    }

    #[test]
    fn polygon_mask_uses_even_odd_rule() {
        let p = rectangular(5, 5, 1.0, 1.0).unwrap();
        let square = [(-1.5, -1.5), (1.5, -1.5), (1.5, 1.5), (-1.5, 1.5)];
        assert_eq!(polygon_mask(&p, &square).unwrap().len(), 9);
        let triangle = [(-3.0, -3.0), (3.5, -3.0), (-3.0, 3.5)];
        let kept = polygon_mask(&p, &triangle).unwrap();
        assert_eq!(kept.len(), 15);
        assert!(kept.iter().all(|q| q.x + q.y <= 1e-12));
        assert!(polygon_mask(&p, &square[..2]).is_err());
        assert!(polygon_mask(&p, &[(0.0, 0.0), (1.0, f64::NAN), (1.0, 1.0)]).is_err());
    }
}
