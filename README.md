# antenna-r

`antenna_r` is a fast, general-purpose Rust library that computes **antenna gain in any direction**.
Every antenna model implements one trait, `GainPattern`, so callers ask for "the gain in dBi towards
this direction" without caring which antenna it is:

- **Parabolic and aperture antennas**: uniform and tapered circular apertures, Gaussian beams, and
  the ITU-R reference patterns S.465-6, S.580-6, S.1428-1, F.699-9, F.1245-3, S.672-4, S.1528-1,
  BO.1213-1, BO.1295-0 and the Radio Regulations Appendix 30/30A (APS30) patterns, with co- and
  cross-polar components where the texts define them.
- **Custom gain from GRD files**: TICRA GRASP `.grd` far-field grids (reader and writer), CSV
  θ/φ tables and 1-D cuts, with bilinear interpolation.
- **Steerable antennas**: any pattern on an az/el, X-over-Y or free mount, with limits, roll,
  pointing error and electronic scan loss.
- **Phased arrays**: element pattern × array factor, lattices, Taylor and Dolph–Chebyshev tapers,
  phase or true-time-delay steering (beam squint), phase quantisation, failures, multi-beam.

It builds on [`rf-r`](https://github.com/Pierre-Penven/rf-r) (dB math, aperture gain) and
[`geometry_r`](https://github.com/Pierre-Penven/geometry-r) (vectors, satellite and ground frames).
The library is pure and thread-safe: no I/O beyond parsing the text you pass in, and no async code.

## Requirements

- Rust 1.88 or later (edition 2024).
- Network access to GitHub at build time: the two sibling crates are git dependencies until they
  are published on crates.io.

## Installation

```toml
[dependencies]
antenna_r = { git = "https://github.com/Pierre-Penven/antenna-r" }
```

Optional feature `parallel` adds rayon-based batch evaluation (`antenna_r::batch`).

## Usage

Conventions: antenna boresight is `+Z`; `θ` is the angle from boresight, `φ` the azimuth from `+X`
towards `+Y`. Units are SI (metres, hertz, radians) and every name carries its unit.

### A dish and its ITU-R envelope

```rust
use antenna_r::prelude::*;

let dish = TaperedAperture::new(1.2, 12e9, 2, -12.0, 0.85).unwrap();
let s465 = ItuS465::from_dish(1.2, 12e9, 0.65, S465Variant::Standard).unwrap();
let ten_degrees = Direction::from_theta_phi(10f64.to_radians(), 0.0);
assert!(dish.gain_dbi(ten_degrees) < s465.gain_dbi(ten_degrees));
assert!((s465.gain_dbi(ten_degrees) - 7.0).abs() < 1e-9); // 32 − 25 log 10
```

### A pattern from a GRD file

```rust
use antenna_r::prelude::*;

let text = "measured beam\n++++\n1\n1 3 2 7\n0 0\n0 0 90 10\n2 2 0\n\
            3.16 0 0 0\n3.16 0 0 0\n1 0 0 0\n1 0 0 0\n";
let grid = GrdFile::parse(text).unwrap().pattern(0).unwrap();
let g = grid.gain_dbi(Direction::from_theta_phi(5f64.to_radians(), 0.0));
assert!((g - 5.0).abs() < 0.01); // halfway between ≈10 dBi and 0 dBi
```

### A steered dish

```rust
use antenna_r::prelude::*;

let beam = GaussianBeam::new(40.0, 1f64.to_radians()).unwrap();
let mut dish = Steered::new(beam, Mount::AzEl);
let target = Direction::from_theta_phi(40f64.to_radians(), 1.0); // mount (ENU) frame
dish.point_at(target).unwrap();
assert!((dish.gain_dbi(target) - 40.0).abs() < 1e-9);
```

### A phased array

```rust
use antenna_r::prelude::*;
use antenna_r::array::{lattice, taper};

let lambda = rf_r::fspl::wavelength_m(20e9);
let positions = lattice::rectangular(16, 16, lambda / 2.0, lambda / 2.0).unwrap();
let w = taper::chebyshev(16, 30.0).unwrap();
let mut array = PhasedArray::new(CosineQ::new(1.0).unwrap(), &positions, 20e9)
    .unwrap()
    .with_amplitudes(&taper::separable(&w, &w))
    .unwrap();
let target = Direction::from_theta_phi(30f64.to_radians(), 0.0);
array.steer(target, SteeringMode::TrueTimeDelay);
assert!(array.gain_dbi(target) > 25.0);
```

More complete programs live in `examples/`:

```bash
cargo run --example dish_pattern
cargo run --example grd_lookup
cargo run --example steered_dish
cargo run --example phased_array_scan
```

## Troubleshooting

- **`AntennaError::Parse { line, .. }` when reading a GRD file**: the record at that line does not
  match the TICRA layout (see `antenna_r::grd::ticra`). Fortran `D` exponents are accepted.
- **`AntennaError::UnsupportedFormat`**: the file uses an `ICOMP` (5–9) or `IGRID` (2, 3, 8+) code
  that is not supported yet.
- **Gains of `-100` dBi**: that is the default floor of patterns with nulls (`DEFAULT_FLOOR_DBI`).
  Change it with `with_floor_dbi`.

## Development

```bash
cargo test --all-features     # unit, integration and doc tests
make check                    # everything CI runs (fmt, clippy, tests, docs, 100 % coverage, audit, MSRV)
cargo bench                   # criterion benchmarks
```

Issues and pull requests: <https://github.com/Pierre-Penven/antenna-r>.

## Documentation

- [Architecture (arc42)](docs/arc42/architecture.md)
- [C4 diagrams](docs/c4/1-context.md)
- [Technical design](docs/technical-design.md)
- [Architecture decision records](docs/adr/)

## Licence

MIT, see [LICENSE](LICENSE).
