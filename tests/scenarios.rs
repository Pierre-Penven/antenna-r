//! End-to-end scenarios across antenna families.

use antenna_r::analysis::{beamwidth_rad, mean_gain_db_piecewise};
use antenna_r::array::{lattice, taper};
use antenna_r::frames::{direction_from_station, gain_from_station};
use antenna_r::grd::{GridAxis, GridType};
use antenna_r::parabolic::itu::BO1295_PEAK_GAIN_DBI;
use antenna_r::prelude::*;
use geometry_r::coords::{Geodetic, LocalFrame};
use geometry_r::satellite::geo_satellite_position;
use std::f64::consts::PI;

fn close(actual: f64, expected: f64, tol: f64) {
    assert!((actual - expected).abs() <= tol, "expected {expected}, got {actual}");
}

#[test]
fn any_pattern_dispatches_to_every_model() {
    let grid = GridPattern::new(
        GridType::Uv,
        GridAxis::new(-0.1, 0.1, 3).unwrap(),
        GridAxis::new(-0.1, 0.1, 3).unwrap(),
        vec![1.0; 9],
    )
    .unwrap();
    let models: Vec<AnyPattern> = vec![
        Isotropic::new(2.0).unwrap().into(),
        CosineQ::new(1.0).unwrap().into(),
        UniformAperture::new(1.0, 12e9, 0.7).unwrap().into(),
        TaperedAperture::new(1.0, 12e9, 1, -10.0, 0.9).unwrap().into(),
        GaussianBeam::new(30.0, 0.05).unwrap().into(),
        ItuS465::new(100.0, 45.0, S465Variant::Standard).unwrap().into(),
        ItuS580::new(100.0, 45.0).unwrap().into(),
        ItuS672::circular(40.0, 0.02, S672SideLobe::Minus20).unwrap().into(),
        ItuS1428::new(50.0).unwrap().into(),
        ItuF699::from_d_over_lambda(150.0, 18e9).unwrap().into(),
        ItuF1245::from_d_over_lambda(150.0, 18e9, F1245Model::Generalized)
            .unwrap()
            .into(),
        ItuS1528::circular(40.0, 0.02, S1528SideLobe::Minus25).unwrap().into(),
        ItuS1528Ngso::circular(35.0, 0.03, S1528Orbit::Leo, 0.0).unwrap().into(),
        ItuS1528Taylor::new(40.0, 20.0, 4, 0.5, 0.4, 12e9).unwrap().into(),
        ItuBo1213::new(23.4, 0.65, PolarComponent::CrossPolar).unwrap().into(),
        ItuBo1295::new(BO1295_PEAK_GAIN_DBI, PolarComponent::CoPolar)
            .unwrap()
            .into(),
        Ap30EarthReceive::new(37.0, 0.035, Ap30EarthKind::Region2, PolarComponent::CoPolar)
            .unwrap()
            .into(),
        Ap30Satellite::circular(
            43.0,
            0.02,
            Ap30SatelliteKind::TransmitRegions1And3,
            PolarComponent::CoPolar,
        )
        .unwrap()
        .into(),
        Ap30aEarthTransmit::new(57.0, Ap30aPlan::Orb88, PolarComponent::CoPolar)
            .unwrap()
            .into(),
        grid.into(),
    ];
    let d = Direction::from_theta_phi(0.01, 0.3);
    for m in &models {
        let direct = match m {
            AnyPattern::Isotropic(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::CosineQ(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::UniformAperture(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::TaperedAperture(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::GaussianBeam(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS465(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS580(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS672(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS1428(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuF699(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuF1245(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS1528(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS1528Ngso(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuS1528Taylor(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuBo1213(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::ItuBo1295(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::Ap30EarthReceive(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::Ap30Satellite(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::Ap30aEarthTransmit(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
            AnyPattern::Grid(p) => (p.gain_dbi(d), p.peak_gain_dbi()),
        };
        assert_eq!((m.gain_dbi(d), m.peak_gain_dbi()), direct);
    }
}

#[test]
fn ku_band_earth_station_sees_adjacent_satellite_below_s465_envelope() {
    // 2.4 m dish at 12 GHz tracking a GSO satellite; neighbour 2° away.
    let station = LocalFrame::new(&Geodetic::from_degrees(48.0, 2.0, 0.0).unwrap());
    let wanted = geo_satellite_position(5f64.to_radians());
    let neighbour = geo_satellite_position(7f64.to_radians());
    let dish = TaperedAperture::new(2.4, 12e9, 2, -12.0, 0.85).unwrap();
    let mut antenna = Steered::new(dish, Mount::AzEl);
    antenna
        .point_at(direction_from_station(&station, wanted).unwrap())
        .unwrap();

    let on_axis = gain_from_station(&antenna, &station, wanted).unwrap();
    close(on_axis, antenna.peak_gain_dbi(), 1e-9);

    let off_axis = gain_from_station(&antenna, &station, neighbour).unwrap();
    let separation_deg = direction_from_station(&station, wanted)
        .unwrap()
        .angle_to(direction_from_station(&station, neighbour).unwrap())
        .to_degrees();
    let envelope = ItuS465::from_dish(2.4, 12e9, 0.65, S465Variant::Standard)
        .unwrap()
        .gain_at_deg(separation_deg);
    assert!(
        off_axis < envelope,
        "{off_axis} dBi above the S.465 envelope {envelope} dBi"
    );
}

#[test]
fn grd_round_trip_preserves_beamwidth_and_efficiency() {
    let dish = UniformAperture::new(1.0, 20e9, 1.0).unwrap();
    let x = GridAxis::new(0.0, 358.0, 180).unwrap();
    let y = GridAxis::new(0.0, 30.0, 3001).unwrap();
    let file = GrdFile::from_pattern(&dish, GridType::ThetaPhi, x, y);
    let parsed = GrdFile::parse(&file.to_grd_string()).unwrap();
    let grid = parsed
        .pattern(0)
        .unwrap()
        .with_out_of_grid(OutOfGrid::Floor(-100.0))
        .unwrap();
    let bw_model = beamwidth_rad(&dish, 0.0, 3.0).unwrap();
    let bw_grid = beamwidth_rad(&grid, 0.0, 3.0).unwrap();
    close(bw_grid, bw_model, 1e-4);
    let edges = [0.0, 0.1, 30f64.to_radians()];
    let model_mean = mean_gain_db_piecewise(&dish, &edges, 200, 8);
    let grid_mean = mean_gain_db_piecewise(&grid, &edges, 200, 8);
    close(grid_mean, model_mean, 0.01);
}

#[test]
fn taylor_planar_array_scans_with_expected_scan_loss() {
    let f = 20e9;
    let lambda = rf_r::fspl::wavelength_m(f);
    let (n, d) = (16, lambda / 2.0);
    let positions = lattice::rectangular(n, n, d, d).unwrap();
    let coords: Vec<f64> = (0..n).map(|i| (i as f64 - 7.5) * d).collect();
    let wx = taper::taylor(&coords, n as f64 * d, 5, 25.0).unwrap();
    let w = taper::separable(&wx, &wx);
    let element = CosineQ::new(1.0).unwrap();
    let mut array = PhasedArray::new(element, &positions, f)
        .unwrap()
        .with_amplitudes(&w)
        .unwrap();
    let broadside = array.gain_dbi(Direction::BORESIGHT);
    let scan = 40f64.to_radians();
    let target = Direction::from_theta_phi(scan, PI / 4.0);
    array.steer(target, SteeringMode::TrueTimeDelay);
    let scanned = array.gain_dbi(target);
    // cos¹ element: the scanned peak drops by 10 log10(cos θ).
    close(scanned - broadside, 10.0 * scan.cos().log10(), 1e-9);
    assert!(!antenna_r::array::grating::has_grating_lobe((d, 0.0), (0.0, d), lambda, target).unwrap());
}

#[test]
fn steered_array_element_in_a_satellite_body_frame() {
    // An electronically steered beam modelled as a single pattern with scan loss.
    let beam = GaussianBeam::new(35.0, 2f64.to_radians()).unwrap();
    let mut s = Steered::new(beam, Mount::Free).with_scan_loss(ScanLoss::new(1.3, true).unwrap());
    let target = Direction::from_theta_phi(45f64.to_radians(), 1.0);
    s.point_at(target).unwrap();
    close(s.gain_dbi(target), 35.0 + 13.0 * 45f64.to_radians().cos().log10(), 1e-9);
}
