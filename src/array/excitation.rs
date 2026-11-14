//! Hardware limits on element excitations: phase-only control, attenuator
//! range and resolution, phase-shifter resolution.
//!
//! An array beamformer sets one amplitude and one phase per element; real
//! hardware restricts both. [`ExcitationLimits`] describes the amplitude side
//! (applied to the amplitudes before building the array); the phase side is
//! the array's own [`PhasedArray::with_phase_bits`](super::PhasedArray::with_phase_bits)
//! (applied to the total phase). [`quantise_phase`] is the same rounding, for
//! callers that work on phases directly (synthesis).

use std::f64::consts::TAU;

use crate::error::{Result, invalid};

/// Amplitude limits of a beamformer.
///
/// Amplitudes are taken relative to the largest one (0 dB). With
/// `phase_only` every element is driven at full amplitude. Otherwise an
/// attenuator covers `range_db` below full amplitude (smaller amplitudes,
/// zero included, are raised to `−range_db`), in steps of
/// `range_db / (2^bits − 1)` when `bits` is set.
///
/// ```
/// use antenna_r::array::excitation::ExcitationLimits;
/// let att = ExcitationLimits { phase_only: false, range_db: Some(30.0), bits: Some(2) };
/// let a = att.apply(&[1.0, 0.5, 0.0]).unwrap();
/// // Levels 0, −10, −20, −30 dB: 0.5 (−6 dB) rounds to −10 dB, 0 is raised to −30 dB.
/// assert!((a[1] - 10f64.powf(-0.5)).abs() < 1e-12);
/// assert!((a[2] - 10f64.powf(-1.5)).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ExcitationLimits {
    /// No amplitude control: every element at full amplitude.
    pub phase_only: bool,
    /// Attenuation range below full amplitude (dB, > 0); `None` = unlimited.
    pub range_db: Option<f64>,
    /// Attenuator resolution (1–16 bits over `range_db`); `None` = continuous.
    pub bits: Option<u8>,
}

impl ExcitationLimits {
    /// Checks the limits.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-positive or non-finite range, `bits` outside 1–16, or `bits`
    /// without a range.
    pub fn validate(&self) -> Result<()> {
        if let Some(r) = self.range_db
            && !(r.is_finite() && r > 0.0)
        {
            return Err(invalid("range_db", format!("must be finite and > 0, got {r}")));
        }
        if let Some(b) = self.bits {
            if !(1..=16).contains(&b) {
                return Err(invalid("bits", format!("must be in 1..=16, got {b}")));
            }
            if self.range_db.is_none() {
                return Err(invalid("bits", "an attenuator resolution needs range_db"));
            }
        }
        Ok(())
    }

    /// Whether the limits change nothing.
    #[must_use]
    pub fn is_unlimited(&self) -> bool {
        !self.phase_only && self.range_db.is_none()
    }

    /// Amplitudes as the hardware can set them, relative to the largest (≤ 1).
    ///
    /// # Errors
    /// As [`validate`](Self::validate); also negative or non-finite
    /// amplitudes, or all zero.
    pub fn apply(&self, amplitudes: &[f64]) -> Result<Vec<f64>> {
        self.validate()?;
        if amplitudes.iter().any(|a| !(a.is_finite() && *a >= 0.0)) {
            return Err(invalid("amplitudes", "must be finite and ≥ 0"));
        }
        let max = amplitudes.iter().copied().fold(0.0, f64::max);
        if max <= 0.0 {
            return Err(invalid("amplitudes", "no radiating element"));
        }
        Ok(amplitudes.iter().map(|a| self.level(a / max)).collect())
    }

    /// One relative amplitude (`0 ≤ a ≤ 1`) as the hardware sets it.
    fn level(&self, a: f64) -> f64 {
        if self.phase_only {
            return 1.0;
        }
        let Some(range) = self.range_db else { return a };
        let db = if a > 0.0 {
            (20.0 * a.log10()).max(-range)
        } else {
            -range
        };
        let db = match self.bits {
            Some(b) => {
                let step = range / f64::from((1u32 << b) - 1);
                ((db / step).round() * step).max(-range)
            }
            None => db,
        };
        10f64.powf(db / 20.0)
    }
}

/// Phase rounded to `bits` bits (`2π / 2^bits` steps); unchanged for `None`.
///
/// ```
/// use antenna_r::array::excitation::quantise_phase;
/// use std::f64::consts::FRAC_PI_2;
/// assert_eq!(quantise_phase(1.0, Some(2)), FRAC_PI_2);
/// assert_eq!(quantise_phase(1.0, None), 1.0);
/// ```
#[must_use]
pub fn quantise_phase(phase_rad: f64, bits: Option<u8>) -> f64 {
    match bits {
        Some(b) => {
            let q = TAU / f64::from(1u32 << b);
            (phase_rad / q).round() * q
        }
        None => phase_rad,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;

    #[test]
    fn unlimited_only_normalises() {
        let l = ExcitationLimits::default();
        assert!(l.is_unlimited());
        assert_eq!(l.apply(&[2.0, 1.0, 0.0]).unwrap(), vec![1.0, 0.5, 0.0]);
    }

    #[test]
    fn phase_only_drives_every_element_fully() {
        let l = ExcitationLimits {
            phase_only: true,
            ..Default::default()
        };
        assert!(!l.is_unlimited());
        assert_eq!(l.apply(&[2.0, 1.0, 0.0]).unwrap(), vec![1.0; 3]);
    }

    #[test]
    fn continuous_attenuator_clamps_to_its_range() {
        let l = ExcitationLimits {
            phase_only: false,
            range_db: Some(20.0),
            bits: None,
        };
        let a = l.apply(&[1.0, 0.5, 0.01, 0.0]).unwrap();
        assert_close(a[1], 0.5, 1e-12);
        assert_close(a[2], 0.1, 1e-12);
        assert_close(a[3], 0.1, 1e-12);
    }

    #[test]
    fn quantised_attenuator_rounds_in_db() {
        let l = ExcitationLimits {
            phase_only: false,
            range_db: Some(3.0),
            bits: Some(1),
        };
        let a = l.apply(&[1.0, 0.8, 0.6]).unwrap();
        // Levels 0 and −3 dB: 0.8 (−1.9 dB) → −3 dB, 0.6 (−4.4 dB) → −3 dB.
        assert_close(a[0], 1.0, 1e-12);
        assert_close(a[1], 10f64.powf(-0.15), 1e-12);
        assert_close(a[2], 10f64.powf(-0.15), 1e-12);
    }

    #[test]
    fn invalid_limits_and_amplitudes_are_rejected() {
        let bad = |l: ExcitationLimits, name: &str| {
            assert!(
                matches!(l.validate(), Err(AntennaError::InvalidParameter { name: n, .. }) if n == name),
                "{l:?}"
            );
        };
        bad(
            ExcitationLimits {
                range_db: Some(0.0),
                ..Default::default()
            },
            "range_db",
        );
        bad(
            ExcitationLimits {
                range_db: Some(f64::NAN),
                ..Default::default()
            },
            "range_db",
        );
        bad(
            ExcitationLimits {
                range_db: Some(10.0),
                bits: Some(0),
                ..Default::default()
            },
            "bits",
        );
        bad(
            ExcitationLimits {
                range_db: Some(10.0),
                bits: Some(17),
                ..Default::default()
            },
            "bits",
        );
        bad(
            ExcitationLimits {
                bits: Some(3),
                ..Default::default()
            },
            "bits",
        );
        let l = ExcitationLimits::default();
        assert!(l.apply(&[1.0, -1.0]).is_err());
        assert!(l.apply(&[1.0, f64::NAN]).is_err());
        assert!(l.apply(&[0.0, 0.0]).is_err());
        assert!(
            ExcitationLimits {
                bits: Some(2),
                ..Default::default()
            }
            .apply(&[1.0])
            .is_err()
        );
    }

    #[test]
    fn phases_round_to_the_resolution() {
        assert_close(quantise_phase(0.8, Some(1)), 0.0, 1e-15);
        assert_close(quantise_phase(2.0, Some(1)), std::f64::consts::PI, 1e-15);
        assert_close(quantise_phase(-0.5, Some(3)), -std::f64::consts::FRAC_PI_4, 1e-15);
    }
}
