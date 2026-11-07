//! Criterion benchmarks of gain evaluation for every antenna family.

use antenna_r::array::lattice::rectangular;
use antenna_r::grd::ticra::GrdFile;
use antenna_r::grd::{GridAxis, GridType};
use antenna_r::prelude::*;
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

const FREQ_HZ: f64 = 12e9;
const N_DIRS: usize = 10_000;

fn directions() -> Vec<Direction> {
    (0..N_DIRS)
        .map(|i| Direction::from_theta_phi((i % 100) as f64 * 1e-3, (i / 100) as f64 * 0.06))
        .collect()
}

fn bench_pattern<P: GainPattern>(c: &mut Criterion, name: &str, pattern: &P) {
    let dirs = directions();
    let mut out = vec![0.0; dirs.len()];
    c.bench_function(name, |b| {
        b.iter(|| pattern.gain_dbi_batch(black_box(&dirs), black_box(&mut out)));
    });
}

fn parabolic(c: &mut Criterion) {
    bench_pattern(
        c,
        "uniform_aperture_10k",
        &UniformAperture::new(1.2, FREQ_HZ, 0.65).unwrap(),
    );
    bench_pattern(
        c,
        "tapered_aperture_10k",
        &TaperedAperture::new(1.2, FREQ_HZ, 2, -10.0, 0.8).unwrap(),
    );
    bench_pattern(
        c,
        "itu_s465_10k",
        &ItuS465::from_dish(1.2, FREQ_HZ, 0.65, S465Variant::Standard).unwrap(),
    );
    bench_pattern(
        c,
        "itu_s672_10k",
        &ItuS672::circular(40.0, 0.01, S672SideLobe::Minus20).unwrap(),
    );
    bench_pattern(c, "itu_f699_10k", &ItuF699::from_d_over_lambda(150.0, 18e9).unwrap());
    bench_pattern(
        c,
        "itu_s1528_taylor_10k",
        &ItuS1528Taylor::new(40.0, 20.0, 4, 0.5, 0.5, FREQ_HZ).unwrap(),
    );
    bench_pattern(
        c,
        "ap30_satellite_fast_roll_off_10k",
        &Ap30Satellite::circular(
            43.0,
            0.02,
            Ap30SatelliteKind::FastRollOffRegions1And3,
            PolarComponent::CoPolar,
        )
        .unwrap(),
    );
}

fn grid(c: &mut Criterion) {
    let dish = UniformAperture::new(1.2, FREQ_HZ, 0.65).unwrap();
    let x = GridAxis::new(-0.2, 0.2, 201).unwrap();
    let file = GrdFile::from_pattern(&dish, GridType::Uv, x, x);
    bench_pattern(c, "grd_uv_bilinear_10k", &file.pattern(0).unwrap());
}

fn steered(c: &mut Criterion) {
    let mut s = Steered::new(UniformAperture::new(1.2, FREQ_HZ, 0.65).unwrap(), Mount::AzEl);
    s.point_at(Direction::from_theta_phi(0.5, 1.0)).unwrap();
    bench_pattern(c, "steered_az_el_10k", &s);
}

fn arrays(c: &mut Criterion) {
    let lambda = rf_r::fspl::wavelength_m(FREQ_HZ);
    let pos = rectangular(32, 32, lambda / 2.0, lambda / 2.0).unwrap();
    let mut a = PhasedArray::new(CosineQ::new(1.2).unwrap(), &pos, FREQ_HZ).unwrap();
    a.steer(Direction::from_theta_phi(0.3, 0.5), SteeringMode::Phase);
    let dirs = directions()[..1000].to_vec();
    let mut out = vec![0.0; dirs.len()];
    c.bench_function("phased_array_1024_elements_1k_dirs", |b| {
        b.iter(|| a.gain_dbi_batch(black_box(&dirs), black_box(&mut out)));
    });
}

criterion_group!(benches, parabolic, grid, steered, arrays);
criterion_main!(benches);
