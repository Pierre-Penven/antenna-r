//! Error type of the crate.
//!
//! Errors are only produced by validating constructors, setters and file
//! parsers. Gain evaluation never returns an error.

use geometry_r::error::GeoError;
use thiserror::Error;

/// Every error the crate can return.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum AntennaError {
    /// A constructor or setter received a value outside its valid domain.
    #[error("invalid parameter `{name}`: {reason}")]
    InvalidParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// Why the value was rejected, including the value.
        reason: String,
    },
    /// A pattern file could not be parsed.
    #[error("parse error at line {line}: {reason}")]
    Parse {
        /// 1-based line number in the input.
        line: usize,
        /// What was expected or found.
        reason: String,
    },
    /// The file is valid but uses an option this crate does not support.
    #[error("unsupported format: {0}")]
    UnsupportedFormat(String),
    /// A pointing or index lies outside the allowed range.
    #[error("out of range: {0}")]
    OutOfRange(String),
    /// A geometry computation from `geometry_r` failed.
    #[error("geometry error: {0}")]
    Geometry(#[from] GeoError),
}

/// Result alias used throughout the crate.
pub type Result<T> = std::result::Result<T, AntennaError>;

/// Builds an [`AntennaError::InvalidParameter`].
pub(crate) fn invalid(name: &'static str, reason: impl Into<String>) -> AntennaError {
    AntennaError::InvalidParameter {
        name,
        reason: reason.into(),
    }
}

/// Requires a finite, strictly positive value.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(invalid(name, format!("must be finite and > 0, got {value}")))
    }
}

/// Requires a finite value.
pub(crate) fn require_finite(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(name, format!("must be finite, got {value}")))
    }
}

/// Requires a finite value or `-inf` (used for floors and levels in dB).
pub(crate) fn require_level_db(name: &'static str, value: f64) -> Result<f64> {
    if value.is_finite() || value == f64::NEG_INFINITY {
        Ok(value)
    } else {
        Err(invalid(name, format!("must be finite or -inf, got {value}")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive_and_rejects_others() {
        assert_eq!(require_positive("x", 2.0), Ok(2.0));
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                require_positive("x", bad),
                Err(AntennaError::InvalidParameter { name: "x", .. })
            ));
        }
    }

    #[test]
    fn require_finite_rejects_nan_and_inf() {
        assert_eq!(require_finite("x", -3.0), Ok(-3.0));
        assert!(require_finite("x", f64::NAN).is_err());
        assert!(require_finite("x", f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn require_level_accepts_minus_infinity_only() {
        assert_eq!(require_level_db("f", f64::NEG_INFINITY), Ok(f64::NEG_INFINITY));
        assert_eq!(require_level_db("f", -40.0), Ok(-40.0));
        assert!(require_level_db("f", f64::INFINITY).is_err());
        assert!(require_level_db("f", f64::NAN).is_err());
    }

    #[test]
    fn errors_display_their_context() {
        assert_eq!(
            invalid("diameter_m", "must be > 0").to_string(),
            "invalid parameter `diameter_m`: must be > 0"
        );
        let parse = AntennaError::Parse {
            line: 7,
            reason: "expected 4 numbers".into(),
        };
        assert_eq!(parse.to_string(), "parse error at line 7: expected 4 numbers");
        assert_eq!(
            AntennaError::UnsupportedFormat("ICOMP 9".into()).to_string(),
            "unsupported format: ICOMP 9"
        );
        assert_eq!(AntennaError::OutOfRange("az".into()).to_string(), "out of range: az");
        let geo: AntennaError = GeoError::NoIntersection.into();
        assert!(geo.to_string().starts_with("geometry error"));
    }
}
