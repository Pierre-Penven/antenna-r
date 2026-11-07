//! # `antenna_r`
//!
//! Fast, general-purpose antenna gain patterns behind one interface,
//! [`GainPattern`]: ask for the gain in dBi towards a [`Direction`] of the
//! antenna frame, whatever the antenna.
//!
//! | Module | Antennas |
//! |---|---|
//! | [`parabolic`] | uniform and tapered circular apertures, Gaussian beams |
//! | [`parabolic::itu`] | ITU-R S.465-6, S.580-6, S.672-4 patterns |
//! | [`grd`] | gridded patterns: TICRA GRASP `.grd` files, CSV tables, 1-D cuts |
//! | [`steer`] | any pattern on a mount (az/el, X/Y, free), limits, scan loss |
//! | [`array`](mod@array) | phased arrays: lattices, tapers, steering, multi-beam |
//! | [`frames`] | gain towards ECEF targets via `geometry_r` frames |
//! | [`analysis`] | cuts, maps, beamwidths, side lobes, directivity integration |
//!
//! Conventions: antenna boresight `+Z`, `θ` off boresight, `φ` from `+X`
//! towards `+Y`; SI units inside (m, Hz, rad), unit suffixes on every name;
//! gains are absolute (dBi).
//!
//! ```
//! use antenna_r::prelude::*;
//! let dish = UniformAperture::new(1.2, 12e9, 0.65).unwrap();
//! let g = dish.gain_dbi(Direction::from_theta_phi(0.5f64.to_radians(), 0.0));
//! assert!(g < dish.peak_gain_dbi());
//! ```

pub mod analysis;
pub mod array;
#[cfg(feature = "parallel")]
pub mod batch;
pub mod direction;
pub mod error;
pub mod frames;
pub mod grd;
pub mod math;
pub mod parabolic;
pub mod pattern;
pub mod steer;

pub use direction::Direction;
pub use error::{AntennaError, Result};
pub use pattern::{AnyPattern, CosineQ, GainPattern, Isotropic};

/// The most used types in one import.
pub mod prelude {
    pub use crate::array::{MultiBeamArray, PhasedArray, SteeringMode};
    pub use crate::direction::Direction;
    pub use crate::error::{AntennaError, Result};
    pub use crate::grd::{GridAxis, GridPattern, GridType, Interpolation, OutOfGrid, ticra::GrdFile};
    pub use crate::parabolic::itu::{
        Ap30EarthKind, Ap30EarthReceive, Ap30Satellite, Ap30SatelliteKind, Ap30aEarthTransmit, Ap30aPlan, F1245Model,
        ItuBo1213, ItuBo1295, ItuF699, ItuF1245, ItuS465, ItuS580, ItuS672, ItuS1428, ItuS1528, ItuS1528Ngso,
        ItuS1528Taylor, PolarComponent, S465Variant, S672SideLobe, S1528Orbit, S1528SideLobe,
    };
    pub use crate::parabolic::{GaussianBeam, TaperedAperture, UniformAperture};
    pub use crate::pattern::{AnyPattern, CosineQ, GainPattern, Isotropic};
    pub use crate::steer::{Mount, ScanLoss, Steered, SteeringLimits};
}

/// Compiles the README code blocks as doctests.
#[doc = include_str!("../README.md")]
#[cfg(doctest)]
pub struct ReadmeDoctests;
