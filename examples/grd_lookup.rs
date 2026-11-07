//! Exports a pattern to a TICRA GRASP grid, reads it back and interpolates it.

use antenna_r::grd::{GridAxis, GridType, Interpolation};
use antenna_r::prelude::*;

fn main() -> antenna_r::Result<()> {
    let dish = UniformAperture::new(0.6, 20e9, 0.7)?;
    let u = GridAxis::new(-0.1, 0.1, 101)?;
    let text = GrdFile::from_pattern(&dish, GridType::Uv, u, u).to_grd_string();
    println!("GRD file: {} lines", text.lines().count());

    let file = GrdFile::parse(&text)?;
    let grid = file
        .pattern(0)?
        .with_interpolation(Interpolation::BilinearDb)
        .with_out_of_grid(OutOfGrid::Floor(-20.0))?;
    for deg in [0.0, 0.3, 1.0, 3.0, 10.0] {
        let d = Direction::from_theta_phi(f64::to_radians(deg), 0.4);
        println!(
            "{deg:5.1} deg: model {:7.2} dBi, grid {:7.2} dBi",
            dish.gain_dbi(d),
            grid.gain_dbi(d)
        );
    }
    Ok(())
}
