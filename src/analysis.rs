//! Pattern analysis for any [`GainPattern`]: cuts, gain maps, beamwidths,
//! nulls, side lobes and power integration over the sphere.

use std::f64::consts::{PI, TAU};

use geometry_r::math::Vec3;

use crate::direction::Direction;
use crate::math::gauss_legendre;
use crate::pattern::GainPattern;

/// Angular step used to bracket a beam edge before bisection (radians).
const BRACKET_STEP_RAD: f64 = 1e-3;
/// Bisection iterations: brackets of 1e-3 rad shrink below 1e-15 rad.
const BISECTION_STEPS: usize = 40;

/// Gain along the cut `φ = phi_rad` at each `θ` in `thetas_rad` (θ may be
/// negative: `(−θ, φ)` is `(θ, φ + π)`).
///
/// # Panics
/// Panics if `thetas_rad.len() != out.len()`.
pub fn cut(pattern: &dyn GainPattern, phi_rad: f64, thetas_rad: &[f64], out: &mut [f64]) {
    assert_eq!(thetas_rad.len(), out.len(), "thetas and out must have the same length");
    for (o, &t) in out.iter_mut().zip(thetas_rad) {
        *o = pattern.gain_dbi(Direction::from_theta_phi(t, phi_rad));
    }
}

/// Gain map on a θ/φ grid, row-major with `φ` varying fastest:
/// `out[i_theta * phis.len() + i_phi]`.
///
/// # Panics
/// Panics if `out.len() != thetas_rad.len() * phis_rad.len()`.
pub fn gain_map_theta_phi(pattern: &dyn GainPattern, thetas_rad: &[f64], phis_rad: &[f64], out: &mut [f64]) {
    assert_eq!(
        out.len(),
        thetas_rad.len() * phis_rad.len(),
        "out must hold thetas × phis values"
    );
    for (row, &t) in out.chunks_exact_mut(phis_rad.len().max(1)).zip(thetas_rad) {
        for (o, &p) in row.iter_mut().zip(phis_rad) {
            *o = pattern.gain_dbi(Direction::from_theta_phi(t, p));
        }
    }
}

/// Gain map on a u/v grid (front hemisphere), `u` varying fastest:
/// `out[i_v * us.len() + i_u]`. Points outside the unit circle are `NaN`.
///
/// # Panics
/// Panics if `out.len() != us.len() * vs.len()`.
pub fn gain_map_uv(pattern: &dyn GainPattern, us: &[f64], vs: &[f64], out: &mut [f64]) {
    assert_eq!(out.len(), us.len() * vs.len(), "out must hold us × vs values");
    for (row, &v) in out.chunks_exact_mut(us.len().max(1)).zip(vs) {
        for (o, &u) in row.iter_mut().zip(us) {
            *o = Direction::from_uv(u, v).map_or(f64::NAN, |d| pattern.gain_dbi(d));
        }
    }
}

/// Off-boresight angle on the half-cut `φ` where the gain first falls
/// `level_db` below the boresight gain, or `None` if it never does.
#[must_use]
pub fn beam_edge_rad(pattern: &dyn GainPattern, phi_rad: f64, level_db: f64) -> Option<f64> {
    let threshold = pattern.gain_dbi(Direction::BORESIGHT) - level_db;
    let below = |t: f64| pattern.gain_dbi(Direction::from_theta_phi(t, phi_rad)) < threshold;
    let mut lo = 0.0;
    let mut hi = BRACKET_STEP_RAD;
    while !below(hi) {
        lo = hi;
        hi += BRACKET_STEP_RAD;
        if hi > PI {
            return None;
        }
    }
    for _ in 0..BISECTION_STEPS {
        let mid = f64::midpoint(lo, hi);
        if below(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(f64::midpoint(lo, hi))
}

/// Full beamwidth (radians) at `level_db` below boresight in the plane
/// `φ = phi_rad`: the sum of the beam edges on both half-cuts (`φ`, `φ + π`).
/// Use `level_db = 3.0` for the half-power beamwidth.
///
/// ```
/// use antenna_r::{analysis::beamwidth_rad, parabolic::GaussianBeam};
/// let p = GaussianBeam::new(30.0, 0.05).unwrap();
/// assert!((beamwidth_rad(&p, 0.0, 3.0).unwrap() - 0.05).abs() < 1e-9);
/// ```
#[must_use]
pub fn beamwidth_rad(pattern: &dyn GainPattern, phi_rad: f64, level_db: f64) -> Option<f64> {
    Some(beam_edge_rad(pattern, phi_rad, level_db)? + beam_edge_rad(pattern, phi_rad + PI, level_db)?)
}

/// Smallest step of the peak refinement (radians).
const PEAK_RESOLUTION_RAD: f64 = 1e-9;

/// Highest gain near `start`: compass search on the sphere (steps along two
/// directions orthogonal to the current one, halved when no step improves)
/// from `step_rad` down to 1e-9 rad. Returns the direction and its gain.
/// Finds the local maximum whose basin contains `start`.
///
/// ```
/// use antenna_r::{Direction, analysis::refine_peak, parabolic::GaussianBeam, steer::{Mount, Steered}};
/// let mut beam = Steered::new(GaussianBeam::new(30.0, 0.05).unwrap(), Mount::Free);
/// beam.point_at(Direction::from_theta_phi(0.2, 0.3)).unwrap();
/// let (d, g) = refine_peak(&beam, Direction::from_theta_phi(0.18, 0.35), 0.01);
/// assert!((d.theta_rad() - 0.2).abs() < 1e-6 && (g - 30.0).abs() < 1e-9);
/// ```
#[must_use]
pub fn refine_peak(pattern: &dyn GainPattern, start: Direction, step_rad: f64) -> (Direction, f64) {
    let mut best = (start, pattern.gain_dbi(start));
    let mut step = step_rad.abs().min(PI);
    while step > PEAK_RESOLUTION_RAD {
        let v = best.0.vec3();
        // Two unit vectors orthogonal to v.
        let a = if v.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let e1 = v.cross(a).normalize();
        let e2 = v.cross(e1);
        let mut moved = false;
        for t in [e1, e1 * -1.0, e2, e2 * -1.0] {
            // A unit vector plus an orthogonal one is never zero.
            let d = Direction::from_vec3(v * step.cos() + t * step.sin()).unwrap_or(best.0);
            let g = pattern.gain_dbi(d);
            if g > best.1 {
                best = (d, g);
                moved = true;
            }
        }
        if !moved {
            step *= 0.5;
        }
    }
    best
}

/// Highest gain over the sphere: the best point of a θ/φ grid of step
/// `step_rad`, refined with [`refine_peak`]. The step must be finer than the
/// narrowest lobe to be sure to find the global maximum. `None` for a
/// non-positive or non-finite step.
///
/// ```
/// use antenna_r::{Direction, analysis::find_peak, parabolic::GaussianBeam, steer::{Mount, Steered}};
/// let mut beam = Steered::new(GaussianBeam::new(30.0, 0.1).unwrap(), Mount::Free);
/// beam.point_at(Direction::from_theta_phi(1.0, 2.0)).unwrap();
/// let (d, g) = find_peak(&beam, 0.02).unwrap();
/// assert!((d.theta_rad() - 1.0).abs() < 1e-6 && (g - 30.0).abs() < 1e-9);
/// ```
#[must_use]
pub fn find_peak(pattern: &dyn GainPattern, step_rad: f64) -> Option<(Direction, f64)> {
    if !(step_rad.is_finite() && step_rad > 0.0) {
        return None;
    }
    let n_theta = (PI / step_rad).ceil() as usize + 1;
    let n_phi = (TAU / step_rad).ceil() as usize;
    let mut best = (Direction::BORESIGHT, pattern.gain_dbi(Direction::BORESIGHT));
    for i in 0..n_theta {
        let theta = (i as f64 * step_rad).min(PI);
        for j in 0..n_phi {
            let d = Direction::from_theta_phi(theta, j as f64 * step_rad);
            let g = pattern.gain_dbi(d);
            if g > best.1 {
                best = (d, g);
            }
        }
    }
    Some(refine_peak(pattern, best.0, step_rad))
}

/// Main-lobe and first side-lobe figures of a half-cut.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CutFeatures {
    /// Angle of the first minimum after the main lobe (radians).
    pub first_null_rad: f64,
    /// Angle of the highest side lobe after the first null (radians).
    pub peak_sidelobe_rad: f64,
    /// Highest side lobe relative to the boresight gain (dB, ≤ 0).
    pub peak_sidelobe_rel_db: f64,
}

/// Finds the first null and the peak side lobe on the half-cut `φ`, sampling
/// `θ ∈ [0, max_theta_rad]` with `samples` points and refining each extremum
/// with a parabola through its neighbours. Returns `None` when the cut has
/// no local minimum (e.g. a monotonic envelope) or `samples < 3`.
#[must_use]
pub fn cut_features(
    pattern: &dyn GainPattern,
    phi_rad: f64,
    max_theta_rad: f64,
    samples: usize,
) -> Option<CutFeatures> {
    if samples < 3 {
        return None;
    }
    let step = max_theta_rad / (samples - 1) as f64;
    let g: Vec<f64> = (0..samples)
        .map(|i| pattern.gain_dbi(Direction::from_theta_phi(i as f64 * step, phi_rad)))
        .collect();
    let null = (1..samples - 1).find(|&i| g[i] <= g[i - 1] && g[i] < g[i + 1])?;
    let lobe = (null + 1..samples - 1).max_by(|&a, &b| g[a].total_cmp(&g[b]))?;
    let (null_offset, _) = parabolic_vertex(g[null - 1], g[null], g[null + 1]);
    let (lobe_offset, lobe_gain) = parabolic_vertex(g[lobe - 1], g[lobe], g[lobe + 1]);
    Some(CutFeatures {
        first_null_rad: (null as f64 + null_offset) * step,
        peak_sidelobe_rad: (lobe as f64 + lobe_offset) * step,
        peak_sidelobe_rel_db: lobe_gain - g[0],
    })
}

/// Vertex `(offset, value)` of the parabola through `(-1, a)`, `(0, b)`, `(1, c)`.
fn parabolic_vertex(a: f64, b: f64, c: f64) -> (f64, f64) {
    let denom = a - 2.0 * b + c;
    if denom == 0.0 || !denom.is_finite() {
        return (0.0, b);
    }
    let offset = 0.5 * (a - c) / denom;
    (offset, b - 0.25 * (a - c) * offset)
}

/// Mean linear gain over the sphere, `∫ G dΩ / 4π`, in dB, integrating `θ`
/// piecewise over the intervals between consecutive `theta_edges_rad` with
/// `nodes_per_interval` Gauss–Legendre nodes (in `cos θ`) and `n_phi`
/// uniformly spaced azimuths.
///
/// For a pattern expressed as gain the result is the radiation efficiency
/// (0 dB when lossless); for directivity data it is 0 dB by construction.
/// Put extra edges around narrow beams (e.g. `[0, 3·θ3dB, π]`).
///
/// # Panics
/// Panics if fewer than two edges are given, or `nodes_per_interval` or
/// `n_phi` is zero.
#[must_use]
pub fn mean_gain_db_piecewise(
    pattern: &dyn GainPattern,
    theta_edges_rad: &[f64],
    nodes_per_interval: usize,
    n_phi: usize,
) -> f64 {
    assert!(theta_edges_rad.len() >= 2, "need at least two theta edges");
    assert!(n_phi > 0, "n_phi must be > 0");
    let (nodes, weights) = gauss_legendre(nodes_per_interval);
    let dphi = TAU / n_phi as f64;
    let mut total = 0.0;
    for edge in theta_edges_rad.windows(2) {
        // Integrate over μ = cos θ between the interval ends.
        let (mu_a, mu_b) = (edge[0].cos(), edge[1].cos());
        let half = 0.5 * (mu_a - mu_b);
        let mid = f64::midpoint(mu_a, mu_b);
        for (x, w) in nodes.iter().zip(&weights) {
            let theta = (mid + half * x).clamp(-1.0, 1.0).acos();
            let ring: f64 = (0..n_phi)
                .map(|k| rf_r::db::from_db(pattern.gain_dbi(Direction::from_theta_phi(theta, k as f64 * dphi))))
                .sum();
            total += w * half * ring * dphi;
        }
    }
    rf_r::db::to_db(total / (4.0 * PI))
}

/// [`mean_gain_db_piecewise`] over the single interval `[0, π]`.
///
/// ```
/// use antenna_r::{analysis::mean_gain_db, Isotropic};
/// let iso = Isotropic::new(0.0).unwrap();
/// assert!(mean_gain_db(&iso, 8, 4).abs() < 1e-12);
/// ```
///
/// # Panics
/// Same as [`mean_gain_db_piecewise`].
#[must_use]
pub fn mean_gain_db(pattern: &dyn GainPattern, n_theta: usize, n_phi: usize) -> f64 {
    mean_gain_db_piecewise(pattern, &[0.0, PI], n_theta, n_phi)
}

/// Directivity (dBi) implied by a pattern: peak gain minus the mean gain,
/// with the same integration parameters as [`mean_gain_db_piecewise`].
///
/// # Panics
/// Same as [`mean_gain_db_piecewise`].
#[must_use]
pub fn directivity_dbi(
    pattern: &dyn GainPattern,
    theta_edges_rad: &[f64],
    nodes_per_interval: usize,
    n_phi: usize,
) -> f64 {
    pattern.peak_gain_dbi() - mean_gain_db_piecewise(pattern, theta_edges_rad, nodes_per_interval, n_phi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::assert_close;
    use crate::pattern::{CosineQ, Isotropic};

    #[test]
    fn peak_search_on_the_sphere() {
        // Lobes peaks at θ = 0; refined from nearby and from the grid.
        let (d, g) = refine_peak(&Lobes, Direction::from_theta_phi(0.2, 1.0), 0.1);
        assert!(d.theta_rad() < 1e-6, "{}", d.theta_rad());
        assert_close(g, 20.0, 1e-6);
        // Near the X axis the other orthogonal basis is used.
        let x = Direction::from_theta_phi(std::f64::consts::FRAC_PI_2, 0.0);
        let (_, gx) = refine_peak(&Isotropic::new(3.0).unwrap(), x, 0.1);
        assert_close(gx, 3.0, 1e-12);
        let (d, g) = find_peak(&Lobes, 0.05).unwrap();
        assert!(d.theta_rad() < 1e-6);
        assert_close(g, 20.0, 1e-6);
        assert!(find_peak(&Lobes, 0.0).is_none());
        assert!(find_peak(&Lobes, f64::NAN).is_none());
        // A peak away from the grid start (boresight) is found from the grid.
        let mut beam = crate::steer::Steered::new(
            crate::parabolic::GaussianBeam::new(30.0, 0.1).unwrap(),
            crate::steer::Mount::Free,
        );
        beam.point_at(Direction::from_theta_phi(1.0, 2.0)).unwrap();
        let (d, g) = find_peak(&beam, 0.03).unwrap();
        assert_close(d.theta_rad(), 1.0, 1e-6);
        assert_close(g, 30.0, 1e-9);
        // A back-hemisphere peak is found too.
        let q = CosineQ::new(1.0).unwrap();
        assert!(find_peak(&q, 0.3).unwrap().0.theta_rad() < 1e-6);
    }

    /// Peaks at θ = 0, null at θ = 0.5, side lobe at θ = 1.0.
    struct Lobes;
    impl GainPattern for Lobes {
        fn gain_dbi(&self, d: Direction) -> f64 {
            let t = d.theta_rad();
            if t < 0.5 {
                20.0 - 80.0 * t
            } else {
                -20.0 - 10.0 * ((t - 1.0) / 0.5).powi(2)
            }
        }
        fn peak_gain_dbi(&self) -> f64 {
            20.0
        }
    }

    #[test]
    fn cut_samples_the_requested_plane() {
        let p = CosineQ::new(1.0).unwrap();
        let mut out = [0.0; 2];
        cut(&p, 0.3, &[0.0, -0.5], &mut out);
        assert_eq!(out[0], p.peak_gain_dbi());
        assert_close(out[1], p.gain_dbi(Direction::from_theta_phi(0.5, 0.3 + PI)), 1e-12);
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn cut_panics_on_length_mismatch() {
        cut(&CosineQ::new(1.0).unwrap(), 0.0, &[0.0], &mut []);
    }

    #[test]
    fn theta_phi_map_is_phi_major() {
        let p = CosineQ::new(1.0).unwrap();
        let mut out = [0.0; 4];
        gain_map_theta_phi(&p, &[0.0, 0.5], &[0.0, 1.0], &mut out);
        assert_eq!(out[0], out[1]);
        assert_close(out[2], p.gain_dbi(Direction::from_theta_phi(0.5, 0.0)), 1e-12);
    }

    #[test]
    #[should_panic(expected = "thetas × phis")]
    fn theta_phi_map_panics_on_bad_buffer() {
        gain_map_theta_phi(&CosineQ::new(1.0).unwrap(), &[0.0], &[0.0, 1.0], &mut [0.0]);
    }

    #[test]
    fn uv_map_marks_invisible_points_nan() {
        let p = Isotropic::new(1.0).unwrap();
        let mut out = [0.0; 2];
        gain_map_uv(&p, &[0.0, 0.9], &[0.9], &mut out);
        assert_eq!(out[0], 1.0);
        assert!(out[1].is_nan());
    }

    #[test]
    #[should_panic(expected = "us × vs")]
    fn uv_map_panics_on_bad_buffer() {
        gain_map_uv(&Isotropic::new(0.0).unwrap(), &[0.0], &[0.0], &mut []);
    }

    #[test]
    fn beam_edge_of_cosine_pattern() {
        // cos θ falls 3.0103 dB at 60°.
        let p = CosineQ::new(1.0).unwrap();
        let edge = beam_edge_rad(&p, 0.0, 10.0 * 2.0_f64.log10()).unwrap();
        assert_close(edge, PI / 3.0, 1e-9);
        assert_close(
            beamwidth_rad(&p, 0.0, 10.0 * 2.0_f64.log10()).unwrap(),
            2.0 * PI / 3.0,
            1e-9,
        );
    }

    #[test]
    fn beam_edge_is_none_for_flat_patterns() {
        let p = Isotropic::new(0.0).unwrap();
        assert_eq!(beam_edge_rad(&p, 0.0, 3.0), None);
        assert_eq!(beamwidth_rad(&p, 0.0, 3.0), None);
    }

    #[test]
    fn cut_features_find_null_and_side_lobe() {
        assert_eq!(Lobes.peak_gain_dbi(), 20.0);
        let f = cut_features(&Lobes, 0.0, 2.0, 2001).unwrap();
        assert_close(f.first_null_rad, 0.5, 2e-3);
        assert_close(f.peak_sidelobe_rad, 1.0, 1e-6);
        assert_close(f.peak_sidelobe_rel_db, -40.0, 1e-6);
    }

    /// Narrow beam towards +X, flat towards −X.
    struct OneSided;
    impl GainPattern for OneSided {
        fn gain_dbi(&self, d: Direction) -> f64 {
            if d.u() > 0.0 { -100.0 * d.theta_rad() } else { 0.0 }
        }
        fn peak_gain_dbi(&self) -> f64 {
            0.0
        }
    }

    #[test]
    fn beamwidth_is_none_when_one_side_never_drops() {
        assert_eq!(OneSided.peak_gain_dbi(), 0.0);
        assert!(beam_edge_rad(&OneSided, 0.0, 3.0).is_some());
        assert_eq!(beamwidth_rad(&OneSided, 0.0, 3.0), None);
    }

    #[test]
    fn cut_features_none_when_the_null_is_the_last_sample() {
        assert_eq!(cut_features(&Lobes, 0.0, 1.0, 3), None);
    }

    #[test]
    fn cut_features_none_without_null_or_samples() {
        let p = CosineQ::new(1.0).unwrap();
        assert_eq!(cut_features(&p, 0.0, 1.0, 100), None);
        assert_eq!(cut_features(&p, 0.0, 1.0, 2), None);
    }

    #[test]
    fn parabolic_vertex_handles_flat_and_curved_points() {
        assert_eq!(parabolic_vertex(1.0, 1.0, 1.0), (0.0, 1.0));
        let (o, v) = parabolic_vertex(0.0, 1.0, 0.0);
        assert_close(o, 0.0, 1e-15);
        assert_close(v, 1.0, 1e-15);
        // y = −(x − 0.25)² sampled at −1, 0, 1.
        let f = |x: f64| -(x - 0.25) * (x - 0.25);
        let (o, v) = parabolic_vertex(f(-1.0), f(0.0), f(1.0));
        assert_close(o, 0.25, 1e-15);
        assert_close(v, 0.0, 1e-15);
        assert_eq!(parabolic_vertex(f64::INFINITY, 0.0, 0.0), (0.0, 0.0));
    }

    #[test]
    fn mean_gain_of_isotropic_and_cosine() {
        assert_close(mean_gain_db(&Isotropic::new(3.0).unwrap(), 4, 3), 3.0, 1e-12);
        let p = CosineQ::new(4.0).unwrap().with_floor_dbi(f64::NEG_INFINITY).unwrap();
        assert_close(mean_gain_db_piecewise(&p, &[0.0, 0.5, PI], 32, 4), 0.0, 1e-6);
        assert_close(directivity_dbi(&p, &[0.0, PI], 64, 4), 10.0 * 10.0_f64.log10(), 1e-4);
    }

    #[test]
    #[should_panic(expected = "two theta edges")]
    fn mean_gain_needs_two_edges() {
        let _ = mean_gain_db_piecewise(&Isotropic::new(0.0).unwrap(), &[0.0], 4, 4);
    }

    #[test]
    #[should_panic(expected = "n_phi")]
    fn mean_gain_needs_azimuths() {
        let _ = mean_gain_db(&Isotropic::new(0.0).unwrap(), 4, 0);
    }
}
