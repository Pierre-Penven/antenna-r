//! 16 x 16 Taylor-tapered array: scan loss, beam squint and phase quantisation.

use antenna_r::array::{lattice, taper};
use antenna_r::prelude::*;

fn main() -> antenna_r::Result<()> {
    let f0 = 20e9;
    let d = rf_r::fspl::wavelength_m(f0) / 2.0;
    let n = 16;
    let positions = lattice::rectangular(n, n, d, d)?;
    let coords: Vec<f64> = (0..n).map(|i| (i as f64 - (n as f64 - 1.0) / 2.0) * d).collect();
    let wx = taper::taylor(&coords, n as f64 * d, 5, 25.0)?;
    let mut array = PhasedArray::new(CosineQ::new(1.2)?, &positions, f0)?
        .with_amplitudes(&taper::separable(&wx, &wx))?
        .with_phase_bits(Some(5))?;
    println!("taper efficiency {:.3}", array.taper_efficiency());
    for scan_deg in [0.0, 20.0, 40.0, 60.0] {
        let target = Direction::from_theta_phi(f64::to_radians(scan_deg), 0.0);
        array.steer(target, SteeringMode::Phase);
        println!(
            "scan {scan_deg:4.1} deg: {:6.2} dBi at f0, {:6.2} dBi at 1.05 f0",
            array.gain_dbi(target),
            array.gain_dbi_at(target, 1.05 * f0)
        );
    }
    Ok(())
}
