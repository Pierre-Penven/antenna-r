# antenna-r technical design

## 1. Overview

`antenna_r` evaluates far-field antenna gain (dBi) towards any direction. All models implement
`GainPattern`; constructors validate and precompute, evaluation is allocation-free `f64` math.
Architecture context: [arc42](arc42/architecture.md), [C4](c4/1-context.md),
decisions in [docs/adr](adr/).

## 2. Components

| Module | Main items | Depends on |
|---|---|---|
| `direction` | `Direction` (`from_theta_phi`, `from_uv`, `from_uv_back`, `from_vec3`, `from_tilts`, `from_el_over_az`, `from_el_and_az`, `from_az_over_el`, inverses, `ScanAngles` interop) | `geometry_r::math`, `geometry_r::satellite` |
| `pattern` | `GainPattern`, `AnyPattern`, `Isotropic`, `CosineQ`, `DEFAULT_FLOOR_DBI` | `rf-r::db` |
| `parabolic` | `UniformAperture`, `TaperedAperture`, `GaussianBeam` | `math`, `rf-r::antenna`, `rf-r::fspl` |
| `parabolic::itu` | `ItuS465` (+ `S465Variant`), `ItuS580`, `ItuS1428`, `ItuF699`, `ItuF1245` (+ `F1245Model`), `ItuS672` (+ `S672SideLobe`), `ItuS1528` (+ `S1528SideLobe`), `ItuS1528Ngso` (+ `S1528Orbit`), `ItuS1528Taylor`, `ItuBo1213`, `ItuBo1295`, `Ap30EarthReceive` (+ `Ap30EarthKind`), `Ap30Satellite` (+ `Ap30SatelliteKind`), `Ap30aEarthTransmit` (+ `Ap30aPlan`), `PolarComponent` | `math`, `rf-r::db` |
| `grd` | `GridPattern`, `GridAxis`, `GridType`, `Interpolation`, `OutOfGrid` | `analysis` |
| `grd::ticra` | `GrdFile`, `GrdSet`, `FieldComponents`, `MAX_GRID_POINTS`, `MAX_SETS` | `grd`, `math` |
| `grd::table` | `parse_theta_phi_table`, `parse_symmetric_cut` | `grd` |
| `steer` | `Steered<P>`, `Mount`, `AxisRange`, `SteeringLimits`, `ScanLoss`, `pointing_loss_db` | `geometry_r::math` |
| `array` | `PhasedArray<E>` (incl. `with_excitation`, `power_efficiency`), `MultiBeamArray<E>`, `SteeringMode` | `math`, `rf-r` |
| `array::excitation` | `ExcitationLimits` (phase-only, attenuator range / bits), `quantise_phase` | — |
| `array::synthesis` | `synthesize` (alternating projections over `PatternBound`s, fixed or maximised coverage level), `place_nulls` (closed-form minimum-change nulls), `coverage_start` (matched start for wide coverages), `SynthesisOptions`, `Synthesis` | `array::excitation`, `math`, `rayon` (feature `parallel`) |
| `array::lattice` | `linear`, `rectangular`, `triangular`, `concentric_rings`, `circular_mask`, `polygon_mask` | — |
| `array::taper` | `uniform`, `cosine_on_pedestal`, `taylor`, `chebyshev`, `separable`, `taper_efficiency`, `power_efficiency` | `math` |
| `array::grating` | `grating_lobes`, `has_grating_lobe` | — |
| `array::kernel` | `ArrayKernel` (`PhasedArray::kernel`): shift-invariant array gain of a planar array, separable 1D tables or a 2D table, phasor-recurrence fill | `math`, `rayon` (feature `parallel`) |
| `frames` | `direction_from_satellite`, `gain_from_satellite`, `gains_from_satellite`, `direction_from_station`, `gain_from_station` | `geometry_r` |
| `analysis` | `cut`, `gain_map_theta_phi`, `gain_map_uv`, `beam_edge_rad`, `beamwidth_rad`, `cut_features`, `mean_gain_db(_piecewise)`, `directivity_dbi`, `find_peak`, `refine_peak` | `math` |
| `math` | `Complex`, `bessel_jn`, `lambda_n`, `gauss_legendre`, `ellipse_radius`, `chebyshev_t` | — |
| `batch` (feature `parallel`) | `par_gain_dbi_batch`, `par_array_gain_dbi_at` | `rayon` |
| `error` | `AntennaError`, `Result` | `thiserror` |

## 3. Data model

No persistence. Key in-memory structures:

- `Direction`: one `Vec3` (unit vector, antenna frame).
- `GridPattern`: grid type, two `GridAxis` (start, end, count), gains in dB in a `Vec<f64>`
  (X fastest; `NaN` = missing), interpolation and out-of-grid policy.
- `GrdFile`: header lines, `FieldComponents`, `GridType`, `Vec<GrdSet>` with two complex
  components per point.
- `PhasedArray<E>`: element pattern, positions as structure of arrays behind an `Arc` (shared by
  `MultiBeamArray` beams), taper, failure mask, amplitude/phase errors, steering, and derived
  amplitudes, projections `r_n·û0` and phases.
- `Steered<P>`: pattern, mount, limits, axes, roll, bias and the combined rotation matrix.

## 4. Interfaces

Rust API only (see `cargo doc`). Text inputs:

**TICRA GRASP grid (`.grd`, ASCII)**: header lines up to `++++`; `KTYPE` (= 1);
`NSET ICOMP NCOMP IGRID`; `NSET` lines `IX IY`; per set `XS YS XE YE`, `NX NY KLIMIT`, then per
row optional `IS IN` and one line per point with `2·NCOMP` reals. Supported `ICOMP` 1–4,
`IGRID` 1, 4, 5, 6, 7, `NCOMP` 2 or 3 (third component ignored). Example:

```text
example beam
++++
1
1 3 2 7
0 0
0.0 0.0 90.0 10.0
2 2 0
3.16 0.0 0.0 0.0
3.16 0.0 0.0 0.0
1.0 0.0 0.0 0.0
1.0 0.0 0.0 0.0
```

**CSV tables**: `θ_deg, φ_deg, gain_dBi` (or `θ_deg, gain_dBi` for a symmetric cut), separated by
commas, semicolons or spaces; `#` comments; optional first header line; any record order, but the
samples must form a complete regular grid.

## 5. Key flows

- **Evaluate a pattern**: `Direction` → model closed form → floor. `NaN` in, `NaN` out.
- **Steered evaluation**: `d_ant = Rᵀ d_mount` → optional scan-plane compression
  (`a ← a·cos θs` along the scan plane) → inner pattern → `+ 10 n log10(cos θs)`.
- **Array evaluation**: `Σ a_n exp(j(k r_n·û + ψ_n − k_TTD r_n·û0))` → `|AF|²/Σa²` → plus element gain → floor.
- **GRD parse** error paths: every record checks its field count and values; counts are bounded
  before allocation; failures return `AntennaError::Parse { line, reason }`.

## 6. Configuration

None at runtime. Cargo feature `parallel` (off by default) enables `batch`.

## 7. Build, test and release

```bash
cargo build
cargo test --all-features          # unit + integration + README/doc tests
make check                         # fmt, clippy -D warnings, tests, rustdoc -D warnings,
                                   # 100 % line coverage, cargo audit, MSRV 1.88, docs checks
cargo bench                        # criterion: benches/antenna_bench.rs
bash scripts/docs/new_adr.sh "t"   # new ADR
```

Release: not published yet. Before a first release, switch the sibling crates to crates.io
versions, update `CHANGELOG.md`, tag, and only run `cargo publish` with the owner's go-ahead.

## 8. Operations

Library: no logging, metrics or health checks. Failure modes are returned as `AntennaError`.

## 9. Security

Pattern files are untrusted input. The GRD parser bounds the total number of points
(`MAX_GRID_POINTS` = 2²⁴) and sets (`MAX_SETS` = 4096), rejects non-finite numbers and
out-of-range counts, and never indexes with unchecked file values. `unsafe` is forbidden.
`cargo audit` runs in CI; Dependabot watches dependencies.

## 10. Known limitations

- ITU-R S.672-4 recommends 2 (shaped beams) is not modelled: it needs the coverage contour.
- ITU-R patterns follow the editions in force in 2026-11; new editions must be re-checked.
- GRD `ICOMP` 5–9 and `IGRID` 2, 3, 8+ are unsupported; only total power is computed (no
  co/cross-polar split yet).
- Analytic apertures are front-hemisphere models (back hemisphere = floor).
- No separable fast path for rectangular arrays; no circular Taylor taper; no mutual coupling.
- Field normalisation of GRD files follows GRASP's 4π W convention (confirmed by the owner).
