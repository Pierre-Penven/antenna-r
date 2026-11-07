//! A ground az/el dish tracking a GSO satellite and its gain towards a neighbour.

use antenna_r::frames::{direction_from_station, gain_from_station};
use antenna_r::prelude::*;
use geometry_r::coords::{Geodetic, LocalFrame};
use geometry_r::satellite::geo_satellite_position;

fn main() -> antenna_r::Result<()> {
    let station = LocalFrame::new(&Geodetic::from_degrees(43.6, 1.44, 150.0)?);
    let wanted = geo_satellite_position(f64::to_radians(13.0));
    // Gaussian main beam; the floor stands in for the side lobes.
    let beam = GaussianBeam::from_dish(1.8, 11.7e9, 0.65, 70.0)?.with_floor_dbi(-10.0)?;
    let mut dish = Steered::new(beam, Mount::AzEl);
    dish.point_at(direction_from_station(&station, wanted)?)?;
    let (az, el) = dish.axes();
    println!("pointing az {:.2} deg, el {:.2} deg", az.to_degrees(), el.to_degrees());
    for lon in [13.0, 14.0, 16.0, 19.0] {
        let sat = geo_satellite_position(f64::to_radians(lon));
        println!(
            "satellite {lon:4.1} E: {:7.2} dBi",
            gain_from_station(&dish, &station, sat)?
        );
    }
    Ok(())
}
