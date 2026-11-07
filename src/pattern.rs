//! The [`GainPattern`] trait shared by every antenna model, the
//! [`AnyPattern`] enum for static dispatch, and two reference patterns.

use std::sync::Arc;

use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite, require_level_db};
use crate::grd::GridPattern;
use crate::parabolic::itu::{
    Ap30EarthReceive, Ap30Satellite, Ap30aEarthTransmit, ItuBo1213, ItuBo1295, ItuF699, ItuF1245, ItuS465, ItuS580,
    ItuS672, ItuS1428, ItuS1528, ItuS1528Ngso, ItuS1528Taylor,
};
use crate::parabolic::{GaussianBeam, TaperedAperture, UniformAperture};

/// Default lowest gain returned by patterns that have nulls, in dBi. Low
/// enough to be negligible in any integral, finite so sums stay finite.
pub const DEFAULT_FLOOR_DBI: f64 = -100.0;

/// Far-field gain of an antenna in any direction of its own frame.
///
/// Implementations are immutable after construction, `Send + Sync`, never
/// allocate in [`gain_dbi`](GainPattern::gain_dbi) and never panic on finite
/// input. Non-finite directions return `NaN`.
pub trait GainPattern: Send + Sync {
    /// Absolute gain in dBi towards `dir` (antenna frame).
    fn gain_dbi(&self, dir: Direction) -> f64;

    /// Maximum gain of the pattern in dBi.
    fn peak_gain_dbi(&self) -> f64;

    /// Gain for many directions into a caller buffer.
    ///
    /// # Panics
    /// Panics if `dirs.len() != out.len()`.
    fn gain_dbi_batch(&self, dirs: &[Direction], out: &mut [f64]) {
        assert_eq!(dirs.len(), out.len(), "dirs and out must have the same length");
        for (o, d) in out.iter_mut().zip(dirs) {
            *o = self.gain_dbi(*d);
        }
    }
}

impl<P: GainPattern + ?Sized> GainPattern for &P {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        (**self).gain_dbi(dir)
    }
    fn peak_gain_dbi(&self) -> f64 {
        (**self).peak_gain_dbi()
    }
    fn gain_dbi_batch(&self, dirs: &[Direction], out: &mut [f64]) {
        (**self).gain_dbi_batch(dirs, out);
    }
}

impl<P: GainPattern + ?Sized> GainPattern for Box<P> {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        (**self).gain_dbi(dir)
    }
    fn peak_gain_dbi(&self) -> f64 {
        (**self).peak_gain_dbi()
    }
    fn gain_dbi_batch(&self, dirs: &[Direction], out: &mut [f64]) {
        (**self).gain_dbi_batch(dirs, out);
    }
}

impl<P: GainPattern + ?Sized> GainPattern for Arc<P> {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        (**self).gain_dbi(dir)
    }
    fn peak_gain_dbi(&self) -> f64 {
        (**self).peak_gain_dbi()
    }
    fn gain_dbi_batch(&self, dirs: &[Direction], out: &mut [f64]) {
        (**self).gain_dbi_batch(dirs, out);
    }
}

/// Constant gain in every direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Isotropic {
    gain_dbi: f64,
}

impl Isotropic {
    /// Isotropic radiator with the given gain (0 dBi for a lossless one).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain.
    pub fn new(gain_dbi: f64) -> Result<Self> {
        Ok(Self {
            gain_dbi: require_finite("gain_dbi", gain_dbi)?,
        })
    }
}

impl GainPattern for Isotropic {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        // Keeps the NaN contract for non-finite directions.
        self.gain_dbi + 0.0 * dir.w()
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.gain_dbi
    }
}

/// `cos^q θ` element pattern in the front hemisphere, a constant back level
/// behind it: `G = G0 + 10 q log10(cos θ)` floored at `floor_dbi`.
///
/// The default peak `G0 = 10 log10(2 (q + 1))` makes the pattern lossless
/// (its directivity with no back radiation).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CosineQ {
    q: f64,
    peak_gain_dbi: f64,
    floor_dbi: f64,
}

impl CosineQ {
    /// Lossless `cos^q θ` pattern (`q ≥ 0`).
    ///
    /// ```
    /// use antenna_r::{CosineQ, Direction, GainPattern};
    /// let e = CosineQ::new(1.0).unwrap();
    /// assert!((e.peak_gain_dbi() - 6.0206).abs() < 1e-4);
    /// assert!((e.gain_dbi(Direction::from_theta_phi(1.0471975511965976, 0.0)) - 3.0103).abs() < 1e-4);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a negative or non-finite `q`.
    pub fn new(q: f64) -> Result<Self> {
        if !(q.is_finite() && q >= 0.0) {
            return Err(invalid("q", format!("must be finite and ≥ 0, got {q}")));
        }
        Ok(Self {
            q,
            peak_gain_dbi: rf_r::db::to_db(2.0 * (q + 1.0)),
            floor_dbi: DEFAULT_FLOOR_DBI,
        })
    }

    /// Overrides the peak gain (e.g. measured element gain).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite gain.
    pub fn with_peak_gain_dbi(mut self, peak_gain_dbi: f64) -> Result<Self> {
        self.peak_gain_dbi = require_finite("peak_gain_dbi", peak_gain_dbi)?;
        Ok(self)
    }

    /// Sets the floor (and back-hemisphere level) in dBi.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// unless finite or `-inf`.
    pub fn with_floor_dbi(mut self, floor_dbi: f64) -> Result<Self> {
        self.floor_dbi = require_level_db("floor_dbi", floor_dbi)?;
        Ok(self)
    }
}

impl GainPattern for CosineQ {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        let w = dir.w();
        if w.is_nan() {
            return f64::NAN;
        }
        if w <= 0.0 {
            return self.floor_dbi;
        }
        // 10 q log10(w), written with ln for speed.
        let g = self.peak_gain_dbi + self.q * rf_r::db::to_db(w);
        g.max(self.floor_dbi)
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

/// Every built-in leaf pattern, for static dispatch without `dyn`.
#[derive(Debug, Clone, PartialEq)]
pub enum AnyPattern {
    /// See [`Isotropic`].
    Isotropic(Isotropic),
    /// See [`CosineQ`].
    CosineQ(CosineQ),
    /// See [`UniformAperture`].
    UniformAperture(UniformAperture),
    /// See [`TaperedAperture`].
    TaperedAperture(TaperedAperture),
    /// See [`GaussianBeam`].
    GaussianBeam(GaussianBeam),
    /// See [`ItuS465`].
    ItuS465(ItuS465),
    /// See [`ItuS580`].
    ItuS580(ItuS580),
    /// See [`ItuS672`].
    ItuS672(ItuS672),
    /// See [`ItuS1428`].
    ItuS1428(ItuS1428),
    /// See [`ItuF699`].
    ItuF699(ItuF699),
    /// See [`ItuF1245`].
    ItuF1245(ItuF1245),
    /// See [`ItuS1528`].
    ItuS1528(ItuS1528),
    /// See [`ItuS1528Ngso`].
    ItuS1528Ngso(ItuS1528Ngso),
    /// See [`ItuS1528Taylor`].
    ItuS1528Taylor(ItuS1528Taylor),
    /// See [`ItuBo1213`].
    ItuBo1213(ItuBo1213),
    /// See [`ItuBo1295`].
    ItuBo1295(ItuBo1295),
    /// See [`Ap30EarthReceive`].
    Ap30EarthReceive(Ap30EarthReceive),
    /// See [`Ap30Satellite`].
    Ap30Satellite(Ap30Satellite),
    /// See [`Ap30aEarthTransmit`].
    Ap30aEarthTransmit(Ap30aEarthTransmit),
    /// See [`GridPattern`].
    Grid(GridPattern),
}

macro_rules! dispatch {
    ($self:ident, $p:ident => $e:expr) => {
        match $self {
            AnyPattern::Isotropic($p) => $e,
            AnyPattern::CosineQ($p) => $e,
            AnyPattern::UniformAperture($p) => $e,
            AnyPattern::TaperedAperture($p) => $e,
            AnyPattern::GaussianBeam($p) => $e,
            AnyPattern::ItuS465($p) => $e,
            AnyPattern::ItuS580($p) => $e,
            AnyPattern::ItuS672($p) => $e,
            AnyPattern::ItuS1428($p) => $e,
            AnyPattern::ItuF699($p) => $e,
            AnyPattern::ItuF1245($p) => $e,
            AnyPattern::ItuS1528($p) => $e,
            AnyPattern::ItuS1528Ngso($p) => $e,
            AnyPattern::ItuS1528Taylor($p) => $e,
            AnyPattern::ItuBo1213($p) => $e,
            AnyPattern::ItuBo1295($p) => $e,
            AnyPattern::Ap30EarthReceive($p) => $e,
            AnyPattern::Ap30Satellite($p) => $e,
            AnyPattern::Ap30aEarthTransmit($p) => $e,
            AnyPattern::Grid($p) => $e,
        }
    };
}

impl GainPattern for AnyPattern {
    #[inline]
    fn gain_dbi(&self, dir: Direction) -> f64 {
        dispatch!(self, p => p.gain_dbi(dir))
    }
    fn peak_gain_dbi(&self) -> f64 {
        dispatch!(self, p => p.peak_gain_dbi())
    }
}

macro_rules! impl_from_for_any {
    ($($ty:ident => $variant:ident),*) => {
        $(impl From<$ty> for AnyPattern {
            fn from(p: $ty) -> Self {
                Self::$variant(p)
            }
        })*
    };
}

impl_from_for_any!(
    Isotropic => Isotropic,
    CosineQ => CosineQ,
    UniformAperture => UniformAperture,
    TaperedAperture => TaperedAperture,
    GaussianBeam => GaussianBeam,
    ItuS465 => ItuS465,
    ItuS580 => ItuS580,
    ItuS672 => ItuS672,
    ItuS1428 => ItuS1428,
    ItuF699 => ItuF699,
    ItuF1245 => ItuF1245,
    ItuS1528 => ItuS1528,
    ItuS1528Ngso => ItuS1528Ngso,
    ItuS1528Taylor => ItuS1528Taylor,
    ItuBo1213 => ItuBo1213,
    ItuBo1295 => ItuBo1295,
    Ap30EarthReceive => Ap30EarthReceive,
    Ap30Satellite => Ap30Satellite,
    Ap30aEarthTransmit => Ap30aEarthTransmit,
    GridPattern => Grid
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::mean_gain_db;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    #[test]
    fn isotropic_is_constant_and_validated() {
        let p = Isotropic::new(2.0).unwrap();
        assert_eq!(p.gain_dbi(Direction::from_theta_phi(2.0, 1.0)), 2.0);
        assert_eq!(p.peak_gain_dbi(), 2.0);
        assert!(p.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
        assert!(matches!(
            Isotropic::new(f64::NAN),
            Err(AntennaError::InvalidParameter { name: "gain_dbi", .. })
        ));
    }

    #[test]
    fn cosine_q_follows_formula_and_floor() {
        let p = CosineQ::new(2.0).unwrap();
        assert_close(p.peak_gain_dbi(), 10.0 * 6.0_f64.log10(), 1e-12);
        let theta: f64 = 0.7;
        let expected = 10.0 * 6.0_f64.log10() + 20.0 * theta.cos().log10();
        assert_close(p.gain_dbi(Direction::from_theta_phi(theta, 0.3)), expected, 1e-12);
        assert_eq!(p.gain_dbi(Direction::from_theta_phi(2.0, 0.0)), DEFAULT_FLOOR_DBI);
        assert_eq!(
            p.gain_dbi(Direction::from_theta_phi(std::f64::consts::FRAC_PI_2 - 1e-9, 0.0)),
            DEFAULT_FLOOR_DBI
        );
        assert!(p.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
    }

    #[test]
    fn cosine_q_is_lossless_by_default() {
        for q in [0.0, 1.0, 3.0] {
            let p = CosineQ::new(q).unwrap().with_floor_dbi(f64::NEG_INFINITY).unwrap();
            assert_close(mean_gain_db(&p, 64, 8), 0.0, 1e-3);
        }
    }

    #[test]
    fn cosine_q_builders_validate() {
        assert!(matches!(
            CosineQ::new(-1.0),
            Err(AntennaError::InvalidParameter { name: "q", .. })
        ));
        assert!(CosineQ::new(f64::INFINITY).is_err());
        let p = CosineQ::new(1.0).unwrap().with_peak_gain_dbi(5.0).unwrap();
        assert_eq!(p.peak_gain_dbi(), 5.0);
        assert!(CosineQ::new(1.0).unwrap().with_peak_gain_dbi(f64::NAN).is_err());
        assert!(CosineQ::new(1.0).unwrap().with_floor_dbi(f64::NAN).is_err());
    }

    #[test]
    fn batch_matches_scalar_and_checks_lengths() {
        let p = CosineQ::new(1.0).unwrap();
        let dirs = [Direction::BORESIGHT, Direction::from_theta_phi(0.5, 0.0)];
        let mut out = [0.0; 2];
        p.gain_dbi_batch(&dirs, &mut out);
        assert_eq!(out[0], p.gain_dbi(dirs[0]));
        assert_eq!(out[1], p.gain_dbi(dirs[1]));
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn batch_panics_on_length_mismatch() {
        let p = Isotropic::new(0.0).unwrap();
        p.gain_dbi_batch(&[Direction::BORESIGHT], &mut []);
    }

    #[test]
    fn smart_pointers_forward_to_the_pattern() {
        let p = CosineQ::new(1.0).unwrap();
        let d = Direction::from_theta_phi(0.4, 0.0);
        let mut out = [0.0];
        let r = &p;
        assert_eq!(GainPattern::gain_dbi(&r, d), p.gain_dbi(d));
        assert_eq!(GainPattern::peak_gain_dbi(&r), p.peak_gain_dbi());
        GainPattern::gain_dbi_batch(&r, &[d], &mut out);
        assert_eq!(out[0], p.gain_dbi(d));
        let b: Box<dyn GainPattern> = Box::new(p);
        assert_eq!(b.gain_dbi(d), p.gain_dbi(d));
        assert_eq!(b.peak_gain_dbi(), p.peak_gain_dbi());
        b.gain_dbi_batch(&[d], &mut out);
        assert_eq!(out[0], p.gain_dbi(d));
        let a: Arc<dyn GainPattern> = Arc::new(p);
        assert_eq!(a.gain_dbi(d), p.gain_dbi(d));
        assert_eq!(a.peak_gain_dbi(), p.peak_gain_dbi());
        a.gain_dbi_batch(&[d], &mut out);
        assert_eq!(out[0], p.gain_dbi(d));
    }
}
