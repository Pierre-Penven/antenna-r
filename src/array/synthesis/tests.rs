use std::f64::consts::PI;

use super::*;
use crate::array::{SteeringMode, lattice};
use crate::math::assert_close;
use crate::pattern::{CosineQ, Isotropic};

const F: f64 = 10e9;

fn lambda() -> f64 {
    rf_r::fspl::wavelength_m(F)
}

fn ula(n: usize) -> PhasedArray<Isotropic> {
    PhasedArray::new(
        Isotropic::new(0.0).unwrap(),
        &lattice::linear(n, lambda() / 2.0).unwrap(),
        F,
    )
    .unwrap()
}

fn dir(theta_deg: f64) -> Direction {
    let t = theta_deg.to_radians();
    if t >= 0.0 {
        Direction::from_theta_phi(t, 0.0)
    } else {
        Direction::from_theta_phi(-t, PI)
    }
}

/// Flat-top sector |θ| ≤ `half` (ripple `ripple` dB), side lobes beyond `edge` at `sll` dB.
fn sector(half: f64, edge: f64, ripple: f64, sll: f64) -> Vec<PatternBound> {
    let mut b = Vec::new();
    let mut t = -half;
    while t <= half + 1e-9 {
        b.push(PatternBound {
            direction: dir(t),
            min: Some(Level::Coverage(0.0)),
            max: Some(Level::Coverage(ripple)),
        });
        t += 1.0;
    }
    let mut t = edge;
    while t <= 89.0 {
        for s in [t, -t] {
            b.push(PatternBound {
                direction: dir(s),
                min: None,
                max: Some(Level::Coverage(sll)),
            });
        }
        t += 1.0;
    }
    b
}

/// Gains of `array` with the synthesised excitation at the bound directions.
fn check_gains<E: GainPattern + Clone>(array: &PhasedArray<E>, s: &Synthesis, bounds: &[PatternBound]) {
    let a = array
        .clone()
        .with_excitation(&s.amplitudes, &s.phases_rad)
        .unwrap()
        .with_floor_dbi(-300.0)
        .unwrap();
    for (b, g) in bounds.iter().zip(&s.gains_dbi) {
        assert_close(a.gain_dbi(b.direction), *g, 1e-6);
    }
}

#[test]
fn side_lobes_are_pushed_below_a_ceiling_at_a_fixed_level() {
    let a = ula(16);
    let mut bounds = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Coverage(0.0)),
        max: None,
    }];
    bounds.extend(sector(0.0, 15.0, 0.0, -25.0).into_iter().filter(|b| b.min.is_none()));
    let opts = SynthesisOptions {
        coverage: CoverageLevel::Fixed(11.0),
        iterations: 400,
        ..Default::default()
    };
    let s = synthesize(&a, &bounds, &opts).unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
    assert_eq!(s.coverage_level_dbi, Some(11.0));
    assert!(s.history_db.len() <= 400 && !s.history_db.is_empty());
    check_gains(&a, &s, &bounds);
    // Uniform start: −13 dB side lobes violate the −25 dB ceiling.
    assert!(s.history_db[0] > 5.0);
    assert_close(s.amplitudes.iter().copied().fold(0.0, f64::max), 1.0, 1e-12);
}

#[test]
fn a_flat_top_sector_is_shaped_and_its_level_maximised() {
    let a = ula(24);
    let bounds = sector(12.0, 22.0, 1.5, -18.0);
    let opts = SynthesisOptions {
        iterations: 600,
        ..Default::default()
    };
    let s = synthesize(&a, &bounds, &opts).unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
    let level = s.coverage_level_dbi.unwrap();
    // A 24° sector cannot reach the pencil-beam gain of 10·log10(24) = 13.8 dBi.
    assert!(level < 13.8 && level > 3.0, "{level}");
    check_gains(&a, &s, &bounds);
    for (b, g) in bounds.iter().zip(&s.gains_dbi) {
        if b.min.is_some() {
            assert!(*g >= level - opts.tolerance_db - 1e-9, "{g} < {level}");
        }
    }
}

#[test]
fn hardware_limits_and_failed_elements_are_honoured() {
    let a = ula(16)
        .with_failed_elements(&[3])
        .unwrap()
        .with_phase_bits(Some(4))
        .unwrap();
    let mut a = a;
    a.steer(dir(10.0), SteeringMode::Phase);
    let bounds = sector(3.0, 15.0, 3.0, -15.0);
    let opts = SynthesisOptions {
        coverage: CoverageLevel::Fixed(8.0),
        iterations: 150,
        limits: ExcitationLimits {
            phase_only: true,
            ..Default::default()
        },
        ..Default::default()
    };
    let s = synthesize(&a, &bounds, &opts).unwrap();
    assert_eq!(s.amplitudes[3], 0.0);
    let q = std::f64::consts::TAU / 16.0;
    let steer = a.steering_phases_rad();
    for (n, (amp, ph)) in s.amplitudes.iter().zip(&s.phases_rad).enumerate() {
        if n != 3 {
            assert_close(*amp, 1.0, 1e-12);
            let total = (ph + steer[n]) / q;
            assert!((total - total.round()).abs() < 1e-6, "element {n}: {total}");
        }
    }
    check_gains(&a, &s, &bounds);
}

#[test]
fn element_patterns_and_absolute_bounds_are_used() {
    let a = PhasedArray::new(
        CosineQ::new(2.0).unwrap(),
        &lattice::linear(8, lambda() / 2.0).unwrap(),
        F,
    )
    .unwrap();
    let bounds = vec![
        PatternBound {
            direction: Direction::BORESIGHT,
            min: Some(Level::Dbi(10.0)),
            max: Some(Level::Dbi(30.0)),
        },
        PatternBound {
            direction: dir(40.0),
            min: None,
            max: Some(Level::Dbi(-10.0)),
        },
    ];
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 200,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(s.coverage_level_dbi, None);
    assert!(s.feasible, "{}", s.max_violation_db);
    check_gains(&a, &s, &bounds);
}

#[test]
fn bounds_met_at_the_start_stop_at_once() {
    let a = ula(4);
    let bounds = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Dbi(0.0)),
        max: None,
    }];
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 50,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(s.history_db, vec![0.0]);
    assert_eq!(s.max_violation_db, 0.0);
    assert_eq!(s.amplitudes, vec![1.0; 4]);
}

#[test]
fn contradictory_bounds_are_reported_infeasible() {
    let a = ula(8);
    let bounds = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Coverage(0.0)),
        max: Some(Level::Coverage(-10.0)),
    }];
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 40,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!s.feasible);
    assert!(s.max_violation_db > 4.0);
    assert!(s.coverage_level_dbi.is_some());
    // Only upper coverage-relative bounds: the start level is the pattern maximum.
    let up = vec![PatternBound {
        direction: dir(30.0),
        min: None,
        max: Some(Level::Coverage(-40.0)),
    }];
    let s = synthesize(
        &a,
        &up,
        &SynthesisOptions {
            iterations: 24,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.coverage_level_dbi.is_some());
}

#[test]
fn invalid_inputs_name_their_parameter() {
    let a = ula(4);
    let ok = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Dbi(0.0)),
        max: None,
    }];
    let name = |bounds: &[PatternBound], opts: SynthesisOptions| match synthesize(&a, bounds, &opts) {
        Err(AntennaError::InvalidParameter { name, .. }) => name,
        other => panic!("{other:?}"),
    };
    let d = SynthesisOptions::default();
    assert_eq!(name(&[], d), "bounds");
    let nan = Direction::from_theta_phi(f64::NAN, 0.0);
    assert_eq!(
        name(
            &[PatternBound {
                direction: nan,
                min: None,
                max: None
            }],
            d
        ),
        "bounds"
    );
    let bad_level = PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Coverage(f64::NAN)),
        max: None,
    };
    assert_eq!(name(&[bad_level], d), "bounds");
    let bad_max = PatternBound {
        direction: Direction::BORESIGHT,
        min: None,
        max: Some(Level::Dbi(f64::INFINITY)),
    };
    assert_eq!(name(&[bad_max], d), "bounds");
    assert_eq!(name(&ok, SynthesisOptions { iterations: 0, ..d }), "iterations");
    assert_eq!(name(&ok, SynthesisOptions { tolerance_db: 0.0, ..d }), "tolerance_db");
    assert_eq!(
        name(
            &ok,
            SynthesisOptions {
                tolerance_db: f64::NAN,
                ..d
            }
        ),
        "tolerance_db"
    );
    assert_eq!(
        name(
            &ok,
            SynthesisOptions {
                coverage: CoverageLevel::Fixed(f64::NAN),
                ..d
            }
        ),
        "coverage"
    );
    let limits = ExcitationLimits {
        bits: Some(3),
        ..Default::default()
    };
    assert_eq!(name(&ok, SynthesisOptions { limits, ..d }), "bits");
}

#[test]
fn nulls_are_placed_with_a_minimum_change() {
    let a = ula(16);
    let nulls = [dir(20.0), dir(-35.0)];
    let (amps, phases) = place_nulls(&a, &nulls).unwrap();
    let b = a
        .clone()
        .with_excitation(&amps, &phases)
        .unwrap()
        .with_floor_dbi(-400.0)
        .unwrap();
    for n in nulls {
        assert!(b.gain_dbi(n) < -100.0, "{}", b.gain_dbi(n));
    }
    // The main beam barely moves.
    assert!((b.gain_dbi(Direction::BORESIGHT) - a.gain_dbi(Direction::BORESIGHT)).abs() < 0.5);
    // Failed elements stay off.
    let f = ula(8).with_failed_elements(&[0]).unwrap();
    let (amps, _) = place_nulls(&f, &[dir(30.0)]).unwrap();
    assert_eq!(amps[0], 0.0);
}

#[test]
fn impossible_nulls_are_rejected() {
    let a = ula(4);
    let name = |a: &PhasedArray<Isotropic>, nulls: &[Direction]| match place_nulls(a, nulls) {
        Err(AntennaError::InvalidParameter { name, reason }) => format!("{name}: {reason}"),
        other => panic!("{other:?}"),
    };
    assert!(name(&a, &[]).contains("at least one"));
    assert!(name(&a, &[Direction::from_theta_phi(f64::NAN, 0.0)]).contains("finite"));
    assert!(name(&a, &[dir(10.0), dir(20.0), dir(30.0), dir(40.0)]).contains("radiating"));
    assert!(name(&a, &[dir(10.0), dir(10.0)]).contains("cannot be met"));
    // A null on the beam of a uniform two-element array removes it.
    assert!(name(&ula(2), &[Direction::BORESIGHT]).contains("whole excitation"));
}

#[test]
fn maximised_level_with_lower_bounds_only_on_the_beam() {
    let a = ula(16);
    let mut bounds = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Coverage(0.0)),
        max: None,
    }];
    bounds.extend(sector(0.0, 15.0, 0.0, -20.0).into_iter().filter(|b| b.min.is_none()));
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 300,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
    // The level is the gain on the beam: close to (below) the uniform 12 dBi.
    let level = s.coverage_level_dbi.unwrap();
    assert_close(level, s.gains_dbi[0], 1e-9);
    assert!(level < 12.05 && level > 10.0, "{level}");
}

#[test]
fn a_planar_array_covers_a_disc_and_protects_a_region() {
    let d = lambda() / 2.0;
    let a = PhasedArray::new(
        Isotropic::new(0.0).unwrap(),
        &lattice::rectangular(12, 12, d, d).unwrap(),
        F,
    )
    .unwrap();
    let mut bounds = Vec::new();
    let step = 0.04;
    for i in -25..=25 {
        for j in -25..=25 {
            let (u, v) = (f64::from(i) * step, f64::from(j) * step);
            let Ok(dir) = Direction::from_uv(u, v) else { continue };
            if u.hypot(v) <= 0.2 {
                bounds.push(PatternBound {
                    direction: dir,
                    min: Some(Level::Coverage(0.0)),
                    max: None,
                });
            } else if (u - 0.6).hypot(v) <= 0.15 {
                bounds.push(PatternBound {
                    direction: dir,
                    min: None,
                    max: Some(Level::Coverage(-25.0)),
                });
            }
        }
    }
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 300,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
    let level = s.coverage_level_dbi.unwrap();
    // Wider than the pencil beam (21.6 dBi), so the edge level is lower.
    assert!(level > 8.0 && level < 21.6, "{level}");
    check_gains(&a, &s, &bounds);
}

#[test]
fn shaping_that_misses_is_retried_at_fixed_levels() {
    // Cosine elements, a covered disc and a protected region 25 dB below:
    // the shaping phase misses, a fixed-level run below it meets the bounds.
    let d = lambda() / 2.0;
    let a = PhasedArray::new(
        CosineQ::new(1.2).unwrap(),
        &lattice::rectangular(12, 12, d, d).unwrap(),
        F,
    )
    .unwrap();
    let mut bounds = Vec::new();
    for i in -25..=25 {
        for j in -25..=25 {
            let (u, v) = (f64::from(i) * 0.05, f64::from(j) * 0.05);
            let Ok(dir) = Direction::from_uv(u, v) else { continue };
            if u.hypot(v) <= 0.2 {
                bounds.push(PatternBound {
                    direction: dir,
                    min: Some(Level::Coverage(0.0)),
                    max: None,
                });
            } else if (u - 0.6).hypot(v) <= 0.15 {
                bounds.push(PatternBound {
                    direction: dir,
                    min: None,
                    max: Some(Level::Coverage(-25.0)),
                });
            }
        }
    }
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 300,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
    check_gains(&a, &s, &bounds);
}

#[test]
fn hopeless_bounds_with_a_tiny_budget_end_in_the_shaping_phase() {
    let a = ula(8);
    let bounds = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Coverage(0.0)),
        max: Some(Level::Coverage(-10.0)),
    }];
    for iterations in [3, 12] {
        let s = synthesize(
            &a,
            &bounds,
            &SynthesisOptions {
                iterations,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(!s.feasible);
    }
}

#[test]
fn a_zero_field_at_a_coverage_sample_is_lifted() {
    // Two elements in antiphase: a null on boresight, where gain is required.
    let a = ula(2).with_excitation(&[1.0, 1.0], &[0.0, PI]).unwrap();
    let bounds = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Dbi(0.0)),
        max: None,
    }];
    let s = synthesize(
        &a,
        &bounds,
        &SynthesisOptions {
            iterations: 20,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
}

#[test]
fn a_coverage_far_wider_than_the_beam_is_filled() {
    // A disc of radius 0.25 in uv (about 29° across) with only a lower bound:
    // the 16 × 16 start is a 6° pencil beam whose weakest coverage sample is
    // tens of dB down. The level must approach the flat-top value, not stay
    // at the pencil beam's weakest sample.
    let d = lambda() / 2.0;
    let a = PhasedArray::new(
        Isotropic::new(0.0).unwrap(),
        &lattice::rectangular(16, 16, d, d).unwrap(),
        F,
    )
    .unwrap();
    let mut bounds = Vec::new();
    for i in -25..=25 {
        for j in -25..=25 {
            let (u, v) = (f64::from(i) * 0.01 * 2.0, f64::from(j) * 0.01 * 2.0);
            if u.hypot(v) <= 0.25 {
                bounds.push(PatternBound {
                    direction: Direction::from_uv(u, v).unwrap(),
                    min: Some(Level::Coverage(0.0)),
                    max: None,
                });
            }
        }
    }
    let opts = SynthesisOptions {
        iterations: 300,
        ..Default::default()
    };
    // From the pencil beam: the power-mean target still lifts the coverage.
    let pencil = synthesize(&a, &bounds, &opts).unwrap();
    let from_pencil = pencil.coverage_level_dbi.unwrap();
    assert!(from_pencil > -10.0, "{from_pencil}");
    // From the matched start: close to the flat top.
    let (amps, phases) = coverage_start(&a, &bounds).unwrap();
    let a = a.with_excitation(&amps, &phases).unwrap();
    let s = synthesize(&a, &bounds, &opts).unwrap();
    assert!(s.feasible, "{}", s.max_violation_db);
    let level = s.coverage_level_dbi.unwrap();
    // Parseval: the gain integrates to (λ/d)² = 4 over a uv cell, so a flat
    // top over the disc (area π·0.25²) is at most 10·log10(4 / 0.196) ≈ 13.1 dBi.
    assert!(level > 10.0 && level < 13.2, "{level}");
    assert!(level > from_pencil + 3.0, "{level} vs {from_pencil}");
    // Shaped well beyond the start: more than the first iterate.
    assert!(s.history_db.len() > 1);
    check_gains(&a, &s, &bounds);
}

#[test]
fn the_coverage_start_needs_a_lower_bound_and_skips_failed_elements() {
    let a = ula(4);
    let ceiling = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: None,
        max: Some(Level::Dbi(0.0)),
    }];
    assert!(coverage_start(&a, &ceiling).is_err());
    assert!(coverage_start(&a, &[]).is_err());
    let cover = vec![PatternBound {
        direction: Direction::BORESIGHT,
        min: Some(Level::Dbi(3.0)),
        max: None,
    }];
    let (amps, phases) = coverage_start(&a, &cover).unwrap();
    // One sample on boresight: uniform, in phase.
    for (am, ph) in amps.iter().zip(&phases) {
        assert_close(*am, 1.0, 1e-9);
        assert_close(*ph, 0.0, 1e-9);
    }
    // Failed elements get no weight.
    let (amps, _) = coverage_start(&ula(4).with_failed_elements(&[1]).unwrap(), &cover).unwrap();
    assert_eq!(amps[1], 0.0);
}
