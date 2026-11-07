//! Principal cut of a 1.2 m Ku-band dish against the ITU-R S.465 envelope.

use antenna_r::analysis::{beamwidth_rad, cut_features};
use antenna_r::prelude::*;

fn main() -> antenna_r::Result<()> {
    let (diameter_m, freq_hz) = (1.2, 12e9);
    let dish = TaperedAperture::new(diameter_m, freq_hz, 2, -12.0, 0.85)?;
    let envelope = ItuS465::from_dish(diameter_m, freq_hz, 0.65, S465Variant::Standard)?;

    println!("peak gain        {:.2} dBi", dish.peak_gain_dbi());
    let bw = beamwidth_rad(&dish, 0.0, 3.0).unwrap_or(f64::NAN);
    println!("-3 dB beamwidth  {:.3} deg", bw.to_degrees());
    if let Some(f) = cut_features(&dish, 0.0, 10f64.to_radians(), 20_001) {
        println!(
            "first null {:.3} deg, peak side lobe {:.2} dB",
            f.first_null_rad.to_degrees(),
            f.peak_sidelobe_rel_db
        );
    }
    println!("\n angle   model   S.465");
    for deg in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 60.0] {
        let d = Direction::from_theta_phi(f64::to_radians(deg), 0.0);
        println!("{deg:6.1} {:7.2} {:7.2}", dish.gain_dbi(d), envelope.gain_dbi(d));
    }
    Ok(())
}
