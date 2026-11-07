//! Amplitude tapers (element excitation magnitudes). Every function returns
//! one non-negative amplitude per element, normalised to a maximum of 1.

use std::f64::consts::{FRAC_PI_2, PI};

use geometry_r::math::Vec3;

use crate::error::{Result, invalid, require_positive};
use crate::math::{Complex, chebyshev_t};

fn normalised(mut a: Vec<f64>) -> Vec<f64> {
    let max = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if max > 0.0 {
        for v in &mut a {
            *v = v.abs() / max;
        }
    }
    a
}

fn sll_ratio(sll_db: f64) -> Result<f64> {
    Ok(rf_r::db::amplitude_from_db(require_positive("sll_db", sll_db)?))
}

/// Uniform excitation of `n` elements.
#[must_use]
pub fn uniform(n: usize) -> Vec<f64> {
    vec![1.0; n]
}

/// Taper efficiency `(Σ a)² / (N Σ a²)` in `(0, 1]` (0 for an empty or
/// all-zero taper).
///
/// ```
/// use antenna_r::array::taper::{taper_efficiency, uniform};
/// assert_eq!(taper_efficiency(&uniform(8)), 1.0);
/// ```
#[must_use]
pub fn taper_efficiency(amplitudes: &[f64]) -> f64 {
    let sum: f64 = amplitudes.iter().sum();
    let power: f64 = amplitudes.iter().map(|a| a * a).sum();
    if power > 0.0 {
        sum * sum / (amplitudes.len() as f64 * power)
    } else {
        0.0
    }
}

/// Radial cosine-on-pedestal taper over the X-Y plane:
/// `a = p + (1 − p) cos^e(π r / (2 R))`, `r` clamped to `R`.
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for a non-positive radius, a pedestal outside `[0, 1]` or a negative or
/// non-finite exponent.
pub fn cosine_on_pedestal(positions: &[Vec3], radius_m: f64, pedestal: f64, exponent: f64) -> Result<Vec<f64>> {
    require_positive("radius_m", radius_m)?;
    if !(0.0..=1.0).contains(&pedestal) {
        return Err(invalid("pedestal", format!("must be in [0, 1], got {pedestal}")));
    }
    if !(exponent.is_finite() && exponent >= 0.0) {
        return Err(invalid("exponent", format!("must be finite and ≥ 0, got {exponent}")));
    }
    Ok(normalised(
        positions
            .iter()
            .map(|p| {
                let r = (p.x.hypot(p.y) / radius_m).min(1.0);
                pedestal + (1.0 - pedestal) * (FRAC_PI_2 * r).cos().powf(exponent)
            })
            .collect(),
    ))
}

/// Taylor line-source taper sampled at element coordinates `coords_m`
/// (e.g. the x of each element) over an aperture of length `length_m`
/// centred on 0:
///
/// `g(x) = 1 + 2 Σ_{m=1}^{n̄−1} F(m) cos(2π m x / L)`, with
/// `F(m) = [(n̄−1)!]² / [(n̄−1+m)! (n̄−1−m)!] · Π_{n=1}^{n̄−1} (1 − m² / (σ² (A² + (n − ½)²)))`,
/// `A = acosh(R)/π`, `R = 10^(SLL/20)`, `σ² = n̄² / (A² + (n̄ − ½)²)`.
///
/// For `N` elements at spacing `d` use `length_m = N d`.
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for `nbar < 2`, a non-positive length or side-lobe level.
pub fn taylor(coords_m: &[f64], length_m: f64, nbar: u32, sll_db: f64) -> Result<Vec<f64>> {
    require_positive("length_m", length_m)?;
    if nbar < 2 {
        return Err(invalid("nbar", format!("must be ≥ 2, got {nbar}")));
    }
    let a = sll_ratio(sll_db)?.acosh() / PI;
    let nb = f64::from(nbar);
    let sigma2 = nb * nb / (a * a + (nb - 0.5).powi(2));
    let coefficients: Vec<f64> = (1..nbar)
        .map(|m| {
            let mf = f64::from(m);
            let mut f = factorial_ratio(nbar - 1, m);
            for n in 1..nbar {
                let zn = sigma2 * (a * a + (f64::from(n) - 0.5).powi(2));
                f *= 1.0 - mf * mf / zn;
            }
            f
        })
        .collect();
    Ok(normalised(
        coords_m
            .iter()
            .map(|x| {
                let s = x / length_m;
                1.0 + 2.0
                    * coefficients
                        .iter()
                        .enumerate()
                        .map(|(i, f)| f * (2.0 * PI * (i + 1) as f64 * s).cos())
                        .sum::<f64>()
            })
            .collect(),
    ))
}

/// `[(k)!]² / ((k + m)! (k − m)!)` computed as a product to avoid overflow.
fn factorial_ratio(k: u32, m: u32) -> f64 {
    // = Π_{i=1}^{m} (k − m + i) / (k + i)
    (1..=m).map(|i| f64::from(k - m + i) / f64::from(k + i)).product()
}

/// Dolph–Chebyshev taper of `n` equally spaced elements with all side lobes
/// `sll_db` below the main beam (exact for spacing ≥ λ/2 at broadside).
///
/// The array factor `T_(n−1)(x0 cos(ψ/2))`, `x0 = cosh(acosh(R)/(n − 1))`,
/// is sampled at `n` points and inverted with a DFT.
///
/// # Errors
/// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
/// for `n < 2` or a non-positive side-lobe level.
pub fn chebyshev(n: usize, sll_db: f64) -> Result<Vec<f64>> {
    if n < 2 {
        return Err(invalid("n", format!("needs at least 2 elements, got {n}")));
    }
    let r = sll_ratio(sll_db)?;
    let order = (n - 1) as u32;
    let x0 = (r.acosh() / f64::from(order)).cosh();
    let nf = n as f64;
    let centre = (nf - 1.0) / 2.0;
    let samples: Vec<f64> = (0..n)
        .map(|k| chebyshev_t(order, x0 * (PI * k as f64 / nf).cos()))
        .collect();
    let weights = (0..n)
        .map(|i| {
            let m = i as f64 - centre;
            let sum = samples.iter().enumerate().fold(Complex::default(), |acc, (k, a)| {
                acc + Complex::from_polar(*a, -m * 2.0 * PI * k as f64 / nf)
            });
            sum.re / nf
        })
        .collect();
    Ok(normalised(weights))
}

/// Product taper `wx[i] · wy[j]` for a rectangular grid, X fastest
/// (matching [`super::lattice::rectangular`]).
#[must_use]
pub fn separable(wx: &[f64], wy: &[f64]) -> Vec<f64> {
    wy.iter().flat_map(|y| wx.iter().map(move |x| x * y)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    /// Peak side lobe (dB) of a linear array factor at λ/2 spacing,
    /// broadside, sampled finely over the visible region.
    fn peak_sidelobe_db(w: &[f64]) -> f64 {
        let n = w.len();
        let af = |u: f64| {
            let s = w.iter().enumerate().fold(Complex::default(), |acc, (i, a)| {
                acc + Complex::from_polar(*a, PI * u * (i as f64 - (n as f64 - 1.0) / 2.0))
            });
            s.norm_sqr()
        };
        let peak = af(0.0);
        let samples: Vec<f64> = (0..=20000).map(|k| af(f64::from(k) / 20000.0)).collect();
        let null = (1..samples.len() - 1)
            .find(|&i| samples[i] < samples[i - 1] && samples[i] <= samples[i + 1])
            .unwrap();
        let max = samples[null..].iter().fold(0.0_f64, |m, v| m.max(*v));
        10.0 * (max / peak).log10()
    }

    #[test]
    fn uniform_and_efficiency() {
        assert_eq!(uniform(3), vec![1.0; 3]);
        assert_close(taper_efficiency(&[1.0, 0.0]), 0.5, 1e-15);
        assert_eq!(taper_efficiency(&[]), 0.0);
        assert_eq!(taper_efficiency(&[0.0, 0.0]), 0.0);
        assert_close(peak_sidelobe_db(&uniform(16)), -13.15, 0.1);
    }

    #[test]
    fn chebyshev_side_lobes_are_equal_to_the_design() {
        for (n, sll) in [(10, 30.0), (11, 25.0), (2, 20.0)] {
            let w = chebyshev(n, sll).unwrap();
            assert_eq!(w.len(), n);
            assert_close(w.iter().copied().fold(0.0, f64::max), 1.0, 1e-12);
            // Symmetric.
            for i in 0..n {
                assert_close(w[i], w[n - 1 - i], 1e-12);
            }
            if n > 2 {
                assert_close(peak_sidelobe_db(&w), -sll, 0.01);
            }
        }
    }

    #[test]
    fn chebyshev_known_weights_for_five_elements() {
        // n = 5, R = 10 (20 dB): weights 1 : 1.6085 : 1.9319 normalised.
        let w = chebyshev(5, 20.0).unwrap();
        assert_close(w[0] / w[2], 1.0 / 1.9319, 1e-3);
        assert_close(w[1] / w[2], 1.6085 / 1.9319, 1e-3);
    }

    #[test]
    fn chebyshev_rejects_invalid_input() {
        assert!(matches!(
            chebyshev(1, 30.0),
            Err(AntennaError::InvalidParameter { name: "n", .. })
        ));
        assert!(matches!(
            chebyshev(4, 0.0),
            Err(AntennaError::InvalidParameter { name: "sll_db", .. })
        ));
    }

    #[test]
    fn taylor_reaches_the_design_side_lobe_level() {
        let n = 64;
        let coords: Vec<f64> = (0..n).map(|i| (f64::from(i) - 31.5) * 0.5).collect();
        let w = taylor(&coords, 32.0, 6, 30.0).unwrap();
        let psl = peak_sidelobe_db(&w);
        assert!(psl < -29.5 && psl > -31.0, "peak side lobe {psl}");
        assert_close(w[31], 1.0, 0.01);
        assert!(w[0] < w[31]);
    }

    #[test]
    fn taylor_rejects_invalid_input() {
        assert!(matches!(
            taylor(&[0.0], 1.0, 1, 30.0),
            Err(AntennaError::InvalidParameter { name: "nbar", .. })
        ));
        assert!(taylor(&[0.0], 0.0, 4, 30.0).is_err());
        assert!(taylor(&[0.0], 1.0, 4, -3.0).is_err());
    }

    #[test]
    fn factorial_ratio_matches_definition() {
        // k = 4, m = 2: (4!)² / (6! 2!) = 576 / 1440 = 0.4
        assert_close(factorial_ratio(4, 2), 0.4, 1e-15);
        assert_close(factorial_ratio(4, 0), 1.0, 1e-15);
    }

    #[test]
    fn cosine_on_pedestal_shapes_the_aperture() {
        let pts = [Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 2.0, 0.0)];
        let a = cosine_on_pedestal(&pts, 1.0, 0.2, 2.0).unwrap();
        assert_close(a[0], 1.0, 1e-15);
        assert_close(a[1], 0.2, 1e-12);
        assert_close(a[2], 0.2, 1e-12);
        assert!(cosine_on_pedestal(&pts, 0.0, 0.2, 2.0).is_err());
        assert!(matches!(
            cosine_on_pedestal(&pts, 1.0, 1.2, 2.0),
            Err(AntennaError::InvalidParameter { name: "pedestal", .. })
        ));
        assert!(matches!(
            cosine_on_pedestal(&pts, 1.0, 0.2, f64::NAN),
            Err(AntennaError::InvalidParameter { name: "exponent", .. })
        ));
        // All-zero taper stays zero (pedestal 0, exponent large at the edge only).
        assert_eq!(normalised(vec![0.0, 0.0]), vec![0.0, 0.0]);
    }

    #[test]
    fn separable_is_x_fastest() {
        assert_eq!(separable(&[1.0, 2.0], &[1.0, 3.0]), vec![1.0, 2.0, 3.0, 6.0]);
    }
}
