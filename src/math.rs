//! Numerical helpers: a minimal complex type, Bessel functions of the first
//! kind, Gauss–Legendre quadrature and small geometric utilities.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};
use std::ops::{Add, AddAssign, Mul, Sub, SubAssign};

/// Minimal complex number (`re + j·im`) used for field components and array
/// factors. Only the operations the crate needs are provided.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// Builds `re + j·im`.
    #[inline]
    #[must_use]
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// `amplitude · exp(j·phase)`.
    #[inline]
    #[must_use]
    pub fn from_polar(amplitude: f64, phase_rad: f64) -> Self {
        let (s, c) = phase_rad.sin_cos();
        Self::new(amplitude * c, amplitude * s)
    }

    /// Squared magnitude `re² + im²` (power of a field component).
    #[inline]
    #[must_use]
    pub fn norm_sqr(self) -> f64 {
        self.re.mul_add(self.re, self.im * self.im)
    }

    /// Phase `atan2(im, re)` in `(-π, π]`.
    #[inline]
    #[must_use]
    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
}

impl Complex {
    /// Zero.
    pub const ZERO: Self = Self::new(0.0, 0.0);

    /// Complex conjugate `re − j·im`.
    ///
    /// ```
    /// use antenna_r::math::Complex;
    /// assert_eq!(Complex::new(1.0, 2.0).conj(), Complex::new(1.0, -2.0));
    /// ```
    #[inline]
    #[must_use]
    pub const fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    /// Product with a real number.
    ///
    /// ```
    /// use antenna_r::math::Complex;
    /// assert_eq!(Complex::new(1.0, -2.0).scale(3.0), Complex::new(3.0, -6.0));
    /// ```
    #[inline]
    #[must_use]
    pub fn scale(self, k: f64) -> Self {
        Self::new(self.re * k, self.im * k)
    }

    /// Magnitude `|z|`.
    ///
    /// ```
    /// use antenna_r::math::Complex;
    /// assert_eq!(Complex::new(3.0, 4.0).abs(), 5.0);
    /// ```
    #[inline]
    #[must_use]
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
}

impl Sub for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
}

impl SubAssign for Complex {
    #[inline]
    fn sub_assign(&mut self, o: Self) {
        self.re -= o.re;
        self.im -= o.im;
    }
}

impl Add for Complex {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}

impl AddAssign for Complex {
    #[inline]
    fn add_assign(&mut self, o: Self) {
        self.re += o.re;
        self.im += o.im;
    }
}

impl Mul for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re.mul_add(o.re, -self.im * o.im),
            self.re.mul_add(o.im, self.im * o.re),
        )
    }
}

/// Argument below which the power series is used for `J_n`; above it the
/// Hankel asymptotic expansion is accurate to better than 1e-12.
const BESSEL_SERIES_LIMIT: f64 = 12.0;
/// Relative size at which series terms are negligible.
const SERIES_EPS: f64 = 1e-17;

/// Bessel function of the first kind `J_n(x)` for small integer orders
/// (`n ≤ 8`), accurate to about 1e-12 absolute for all finite `x`.
///
/// Uses the power series `Σ (−1)^k (x/2)^(2k+n) / (k! (k+n)!)` for
/// `|x| < 12` and the Hankel asymptotic expansion otherwise
/// (Abramowitz & Stegun 9.1.10 and 9.2.5).
///
/// ```
/// use antenna_r::math::bessel_jn;
/// assert!((bessel_jn(1, 1.0) - 0.440_050_585_744_933_5).abs() < 1e-12);
/// ```
#[must_use]
pub fn bessel_jn(n: u32, x: f64) -> f64 {
    let ax = x.abs();
    let value = if ax < BESSEL_SERIES_LIMIT {
        bessel_series(n, ax)
    } else {
        bessel_asymptotic(n, ax)
    };
    if x < 0.0 && n % 2 == 1 { -value } else { value }
}

fn bessel_series(n: u32, x: f64) -> f64 {
    let half = 0.5 * x;
    let mut term = 1.0;
    for k in 1..=n {
        term *= half / f64::from(k);
    }
    let q = -half * half;
    let mut sum = term;
    let mut k = 1.0;
    while term.abs() > SERIES_EPS * sum.abs().max(SERIES_EPS) {
        term *= q / (k * (k + f64::from(n)));
        sum += term;
        k += 1.0;
    }
    sum
}

fn bessel_asymptotic(n: u32, x: f64) -> f64 {
    let mu = 4.0 * f64::from(n * n);
    let eight_x = 8.0 * x;
    // t_k = Π_{i=1..k} (μ − (2i−1)²) / (k! (8x)^k); P takes even k with
    // alternating signs, Q the odd ones.
    let mut p = 1.0;
    let mut q = 0.0;
    let mut term = 1.0_f64;
    let mut previous = f64::INFINITY;
    for k in 1..60 {
        let odd = f64::from(2 * k - 1);
        term *= (mu - odd * odd) / (f64::from(k) * eight_x);
        if term.abs() >= previous || term == 0.0 {
            break;
        }
        previous = term.abs();
        let sign = if (k / 2) % 2 == 0 { 1.0 } else { -1.0 };
        if k % 2 == 0 {
            p += sign * term;
        } else {
            q += sign * term;
        }
    }
    let chi = x - (f64::from(n) * FRAC_PI_2 + FRAC_PI_4);
    let (s, c) = chi.sin_cos();
    (2.0 / (PI * x)).sqrt() * (p * c - q * s)
}

/// `Λ_n(x) = 2^(n+1) (n+1)! J_(n+1)(x) / x^(n+1)`, the normalised pattern of
/// the circular aperture distribution `(1 − ρ²)^n`, with `Λ_n(0) = 1`.
/// `Λ_0(x) = 2 J_1(x)/x` is the uniform-aperture ("jinc") pattern.
///
/// ```
/// use antenna_r::math::lambda_n;
/// assert_eq!(lambda_n(0, 0.0), 1.0);
/// assert!(lambda_n(0, 3.831_705_970_207_512).abs() < 1e-12); // first null
/// ```
#[must_use]
pub fn lambda_n(n: u32, x: f64) -> f64 {
    let ax = x.abs();
    if ax < BESSEL_SERIES_LIMIT {
        // Series of J_(n+1)(x)/x^(n+1) scaled: Σ (−1)^k (x/2)^(2k) (n+1)! / (k! (k+n+1)!).
        let m = f64::from(n + 1);
        let q = -0.25 * ax * ax;
        let mut term = 1.0_f64;
        let mut sum = 1.0_f64;
        let mut k = 1.0;
        while term.abs() > SERIES_EPS * sum.abs().max(SERIES_EPS) {
            term *= q / (k * (k + m));
            sum += term;
            k += 1.0;
        }
        sum
    } else {
        let mut scale = 1.0;
        for k in 1..=n + 1 {
            scale *= 2.0 * f64::from(k) / ax;
        }
        scale * bessel_jn(n + 1, ax)
    }
}

/// Gauss–Legendre nodes and weights on `[-1, 1]` (Newton iteration on the
/// Legendre polynomial, accurate to machine precision).
///
/// # Panics
/// Panics if `n == 0`.
#[must_use]
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    assert!(n > 0, "gauss_legendre needs at least one node");
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    let nf = n as f64;
    for i in 0..n.div_ceil(2) {
        let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        for _ in 0..100 {
            let (p, d) = legendre_with_derivative(n, x);
            let dx = p / d;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        let (_, d) = legendre_with_derivative(n, x);
        let w = 2.0 / ((1.0 - x * x) * d * d);
        nodes[i] = -x;
        nodes[n - 1 - i] = x;
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    (nodes, weights)
}

fn legendre_with_derivative(n: usize, x: f64) -> (f64, f64) {
    let mut p0 = 1.0;
    let mut p1 = x;
    for k in 2..=n {
        let kf = k as f64;
        let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
        p0 = p1;
        p1 = p2;
    }
    let pn = if n == 0 { 1.0 } else { p1 };
    let pm1 = if n == 0 { 0.0 } else { p0 };
    let d = n as f64 * (x * pn - pm1) / (x * x - 1.0);
    (pn, d)
}

/// Radius of an ellipse with semi-axes `a` (along angle 0) and `b` at polar
/// angle `angle_rad`: `a b / √((b cos α)² + (a sin α)²)`.
#[inline]
#[must_use]
pub fn ellipse_radius(a: f64, b: f64, angle_rad: f64) -> f64 {
    let (s, c) = angle_rad.sin_cos();
    a * b / (b * c).hypot(a * s)
}

/// Chebyshev polynomial of the first kind `T_n(x)` for any real `x`
/// (`cos(n acos x)` inside `[-1, 1]`, `±cosh(n acosh |x|)` outside).
#[must_use]
pub fn chebyshev_t(n: u32, x: f64) -> f64 {
    let nf = f64::from(n);
    if x.abs() <= 1.0 {
        (nf * x.acos()).cos()
    } else {
        let v = (nf * x.abs().acosh()).cosh();
        if x < 0.0 && n % 2 == 1 { -v } else { v }
    }
}

#[cfg(test)]
pub(crate) fn assert_close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "expected {expected}, got {actual} (tol {tol})"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_subtraction_conjugate_and_scaling() {
        let mut z = Complex::new(3.0, 4.0) - Complex::new(1.0, 1.0);
        assert_eq!(z, Complex::new(2.0, 3.0));
        z -= Complex::new(2.0, 3.0);
        assert_eq!(z, Complex::ZERO);
        assert_eq!(Complex::new(1.0, 2.0).conj().scale(2.0), Complex::new(2.0, -4.0));
        assert_eq!(Complex::new(-3.0, 4.0).abs(), 5.0);
    }

    #[test]
    fn complex_arithmetic_matches_definitions() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, -1.0);
        assert_eq!(a + b, Complex::new(4.0, 1.0));
        assert_eq!(a * b, Complex::new(5.0, 5.0));
        let mut c = a;
        c += b;
        assert_eq!(c, Complex::new(4.0, 1.0));
        assert_close(a.norm_sqr(), 5.0, 1e-15);
        let p = Complex::from_polar(2.0, FRAC_PI_2);
        assert_close(p.re, 0.0, 1e-15);
        assert_close(p.im, 2.0, 1e-15);
        assert_close(p.arg(), FRAC_PI_2, 1e-15);
        assert_eq!(Complex::default(), Complex::new(0.0, 0.0));
    }

    #[test]
    fn bessel_matches_tabulated_values() {
        // Abramowitz & Stegun table 9.1 values.
        assert_close(bessel_jn(0, 0.0), 1.0, 1e-15);
        assert_close(bessel_jn(1, 0.0), 0.0, 1e-15);
        assert_close(bessel_jn(0, 1.0), 0.765_197_686_557_966_6, 1e-13);
        assert_close(bessel_jn(1, 1.0), 0.440_050_585_744_933_5, 1e-13);
        assert_close(bessel_jn(0, 10.0), -0.245_935_764_451_348_3, 1e-12);
        assert_close(bessel_jn(1, 10.0), 0.043_472_746_168_861_44, 1e-12);
    }

    #[test]
    fn bessel_zeros_are_found() {
        assert_close(bessel_jn(0, 2.404_825_557_695_773), 0.0, 1e-12);
        assert_close(bessel_jn(1, 3.831_705_970_207_512), 0.0, 1e-12);
        assert_close(bessel_jn(1, 7.015_586_669_815_619), 0.0, 1e-12);
        assert_close(bessel_jn(1, 13.323_691_936_314_22), 0.0, 1e-12);
        assert_close(bessel_jn(1, 16.470_630_050_877_63), 0.0, 1e-12);
    }

    #[test]
    fn series_and_asymptotic_agree_at_the_switch_point() {
        for n in 0..=4 {
            for x in [12.0, 13.0] {
                let s = bessel_series(n, x);
                let a = bessel_asymptotic(n, x);
                assert_close(s, a, 1e-10);
            }
        }
    }

    #[test]
    fn bessel_satisfies_recurrence_for_large_arguments() {
        // J_(n-1) + J_(n+1) = (2n/x) J_n
        for x in [0.5, 5.0, 25.0, 300.0, 5000.0] {
            for n in 1..=3 {
                let lhs = bessel_jn(n - 1, x) + bessel_jn(n + 1, x);
                let rhs = 2.0 * f64::from(n) / x * bessel_jn(n, x);
                assert_close(lhs, rhs, 1e-11);
            }
        }
    }

    #[test]
    fn bessel_parity_for_negative_arguments() {
        assert_close(bessel_jn(1, -1.0), -bessel_jn(1, 1.0), 1e-15);
        assert_close(bessel_jn(0, -1.0), bessel_jn(0, 1.0), 1e-15);
        assert_close(bessel_jn(2, -20.0), bessel_jn(2, 20.0), 1e-15);
    }

    #[test]
    fn bessel_propagates_nan() {
        assert!(bessel_jn(1, f64::NAN).is_nan());
    }

    #[test]
    fn lambda_matches_its_bessel_definition() {
        for n in 0..=2 {
            for x in [0.3, 2.0, 11.9, 12.1, 40.0, -40.0] {
                let mut scale = 1.0;
                for k in 1..=n + 1 {
                    scale *= 2.0 * f64::from(k) / x;
                }
                assert_close(lambda_n(n, x), scale * bessel_jn(n + 1, x), 1e-11);
            }
            assert_close(lambda_n(n, 0.0), 1.0, 1e-15);
        }
    }

    #[test]
    fn gauss_legendre_integrates_polynomials_exactly() {
        let (x, w) = gauss_legendre(5);
        let integral = |f: &dyn Fn(f64) -> f64| x.iter().zip(&w).map(|(xi, wi)| wi * f(*xi)).sum::<f64>();
        assert_close(integral(&|_| 1.0), 2.0, 1e-14);
        assert_close(integral(&|t| t.powi(8)), 2.0 / 9.0, 1e-14);
        assert_close(integral(&|t| t.powi(3)), 0.0, 1e-14);
    }

    #[test]
    fn gauss_legendre_handles_one_node_and_odd_counts() {
        let (x, w) = gauss_legendre(1);
        assert_close(x[0], 0.0, 1e-15);
        assert_close(w[0], 2.0, 1e-15);
        let (x, w) = gauss_legendre(64);
        assert_close(w.iter().sum::<f64>(), 2.0, 1e-13);
        assert!(x.windows(2).all(|p| p[0] < p[1]));
    }

    #[test]
    #[should_panic(expected = "at least one node")]
    fn gauss_legendre_panics_on_zero_nodes() {
        let _ = gauss_legendre(0);
    }

    #[test]
    fn legendre_of_order_zero_is_one() {
        assert_eq!(legendre_with_derivative(0, 0.3).0, 1.0);
    }

    #[test]
    fn ellipse_radius_on_axes_and_circle() {
        assert_close(ellipse_radius(2.0, 1.0, 0.0), 2.0, 1e-15);
        assert_close(ellipse_radius(2.0, 1.0, FRAC_PI_2), 1.0, 1e-15);
        assert_close(ellipse_radius(1.5, 1.5, 0.7), 1.5, 1e-15);
    }

    #[test]
    fn chebyshev_inside_and_outside_unit_interval() {
        assert_close(chebyshev_t(3, 0.5), 4.0 * 0.125 - 3.0 * 0.5, 1e-15);
        assert_close(chebyshev_t(3, 2.0), 4.0 * 8.0 - 6.0, 1e-12);
        assert_close(chebyshev_t(3, -2.0), -26.0, 1e-12);
        assert_close(chebyshev_t(2, -2.0), 7.0, 1e-12);
    }
}
