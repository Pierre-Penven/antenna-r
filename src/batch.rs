//! Parallel batch evaluation (feature `parallel`, rayon). Results are
//! identical to the serial [`GainPattern::gain_dbi_batch`].

use rayon::prelude::*;

use crate::array::PhasedArray;
use crate::direction::Direction;
use crate::pattern::GainPattern;

/// Directions per rayon task: large enough to amortise scheduling.
const CHUNK: usize = 1024;

/// Gains for many directions on all cores.
///
/// # Panics
/// Panics if `dirs.len() != out.len()`.
pub fn par_gain_dbi_batch<P: GainPattern + ?Sized>(pattern: &P, dirs: &[Direction], out: &mut [f64]) {
    assert_eq!(dirs.len(), out.len(), "dirs and out must have the same length");
    out.par_chunks_mut(CHUNK)
        .zip(dirs.par_chunks(CHUNK))
        .for_each(|(o, d)| pattern.gain_dbi_batch(d, o));
}

/// Phased-array gains at `freq_hz` for many directions on all cores.
///
/// # Panics
/// Panics if `dirs.len() != out.len()`.
pub fn par_array_gain_dbi_at<E: GainPattern>(
    array: &PhasedArray<E>,
    dirs: &[Direction],
    freq_hz: f64,
    out: &mut [f64],
) {
    assert_eq!(dirs.len(), out.len(), "dirs and out must have the same length");
    out.par_iter_mut()
        .zip(dirs.par_iter())
        .with_min_len(CHUNK / 16)
        .for_each(|(o, d)| *o = array.gain_dbi_at(*d, freq_hz));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::array::lattice::rectangular;
    use crate::pattern::{CosineQ, Isotropic};

    fn dirs(n: usize) -> Vec<Direction> {
        (0..n)
            .map(|i| Direction::from_theta_phi(i as f64 * 1e-3, i as f64 * 7e-3))
            .collect()
    }

    #[test]
    fn parallel_batch_equals_serial() {
        let p = CosineQ::new(2.0).unwrap();
        let d = dirs(5000);
        let mut serial = vec![0.0; d.len()];
        let mut parallel = vec![0.0; d.len()];
        p.gain_dbi_batch(&d, &mut serial);
        par_gain_dbi_batch(&p, &d, &mut parallel);
        assert_eq!(serial, parallel);
    }

    #[test]
    fn parallel_array_equals_serial() {
        let lambda = rf_r::fspl::wavelength_m(10e9);
        let pos = rectangular(8, 8, lambda / 2.0, lambda / 2.0).unwrap();
        let a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &pos, 10e9).unwrap();
        let d = dirs(3000);
        let mut out = vec![0.0; d.len()];
        par_array_gain_dbi_at(&a, &d, 11e9, &mut out);
        for (o, dir) in out.iter().zip(&d) {
            assert_eq!(*o, a.gain_dbi_at(*dir, 11e9));
        }
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn parallel_batch_panics_on_length_mismatch() {
        par_gain_dbi_batch(&Isotropic::new(0.0).unwrap(), &dirs(2), &mut [0.0]);
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn parallel_array_panics_on_length_mismatch() {
        let a = PhasedArray::new(Isotropic::new(0.0).unwrap(), &[geometry_r::math::Vec3::ZERO], 1e9).unwrap();
        par_array_gain_dbi_at(&a, &dirs(2), 1e9, &mut [0.0]);
    }
}
