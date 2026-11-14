# antenna-r — CLAUDE.md

Guidance for Claude Code in this repository. It is also the **specification**
of the crate. Milestones M0–M6 are implemented (see §8 for status and the
deliberate deviations); keep this file in sync when the API or scope changes.

`antenna_r` is a fast, general-purpose Rust library that computes **antenna gain
in any direction** for four antenna families behind one common interface:

1. **Parabolic / aperture** antennas (analytic patterns and ITU-R envelopes)
2. **Custom gain** read from **GRD** pattern files (gridded measured or simulated data)
3. **Steerable** antennas (any pattern, pointed mechanically or electronically)
4. **Phased arrays** (element pattern × array factor, with tapers and steering)

The library has no network I/O and no async code. It is pure, deterministic
and thread-safe (file parsing takes `&str` / `impl Read`). It is general-purpose:
do not mention a specific application or customer in code, docs or metadata.

---

## 1. Goals and non-goals

### Goals
- **One interface for every antenna.** Every model implements the
  `GainPattern` trait (§4.1). Callers ask "gain in dBi towards this direction"
  without knowing the antenna type. Models compose: a steerable antenna wraps
  any pattern; a phased array uses any pattern as its element.
- **High performance.** Patterns are built (validated, precomputed) once and
  evaluated millions of times. Evaluation does not allocate, does not return
  `Result` and has a batch variant over slices. Optional `rayon` parallelism.
- **Physical correctness.** Gain is absolute (dBi), normalised consistently
  across models. Every formula is documented with its source and checked
  against reference values in tests.
- **Explicit units.** SI inside (m, Hz, rad, W). Degrees and GHz only in
  explicitly suffixed convenience constructors. Names carry units
  (`_dbi`, `_db`, `_rad`, `_deg`, `_hz`, `_m`).
- **Reuse, don't duplicate.** dB math and basic aperture formulas come from
  `rf-r`; vectors, rotations and scan-angle conventions come from `geometry_r`
  (§2.1).

### Non-goals (for now)
- Full-wave electromagnetic simulation, mutual coupling, active element
  patterns computed from geometry (an *active element pattern* may be supplied
  as a GRD element pattern instead).
- Near-field patterns. Everything is far field.
- Rendering or plotting. The crate exports numbers (cuts, grids); plotting is
  the caller's job.
- Link budgets, path loss, geometry of satellites and ground stations: those
  live in `rf-r` and `geometry_r`.

---

## 2. Toolchain, dependencies and constraints

- **Edition 2024, MSRV Rust 1.88** (`rust-version = "1.88"`). 1.88 is imposed
  by `geometry_r`. Raising it needs an ADR.
- Crate name `antenna_r` (like `geometry_r`; ask before publishing if the owner
  prefers `antenna-r`). Licence MIT, like the sibling crates.
- Not published to crates.io yet. Never run `cargo publish` without the
  owner's explicit go-ahead.

### 2.1 Sibling crates (owner: Pierre-Penven)

Neither is on crates.io yet, so use git dependencies (pin a `rev` or `tag`);
switch to versions once they are published.

```toml
[dependencies]
rf-r       = { git = "https://github.com/Pierre-Penven/rf-r", rev = "d289d8c…" }
geometry_r = { git = "https://github.com/Pierre-Penven/geometry-r", rev = "051bb47…" }
```

Note the real rf-r names: `fspl::wavelength_m`, `constants::SPEED_OF_LIGHT`.

Local checkouts live next to this repo (`../rf-r`, `../geometry-r`). For
co-development use a `[patch]` section in a local, uncommitted config, never a
committed `path` dependency.

**Reuse from `rf-r`** (zero-dependency, `f64` in / `f64` out):
- `db::{to_db, from_db, amplitude_to_db, amplitude_from_db, log_sum, log_sum_all, log_diff}`
  and `*_in_place` batch variants. Do not re-implement dB math.
- `antenna::{parabolic_gain_dbi, aperture_gain_dbi, half_power_beamwidth_deg, g_over_t_db_k}`
  for peak gain and the 70·λ/D rule of thumb.
- `fspl::wavelength_m` and `constants::SPEED_OF_LIGHT`.

**Reuse from `geometry_r`**:
- `math::{Vec3, Mat3}` (rotations `rot_x/y/z`, `mul_vec`, `mul_vec_t`), `wrap_pi`, `wrap_two_pi`.
- `satellite::ScanAngles` — the θ/φ, u/v and X/Y-tilt conventions. Our
  `Direction` type (§4.1) uses exactly the same conventions and converts
  to/from `ScanAngles` losslessly.
- `satellite::SatelliteFrame` and `coords::LocalFrame` — only in the
  integration layer (`frames` module, §4.6), to turn an ECEF target into an
  antenna-frame direction.
- `geometry_r::footprint::gaussian_gain_db` is a simple main-lobe helper;
  full patterns (including the ITU-R envelopes listed in geometry-r's roadmap)
  belong **here**. Do not add pattern code to geometry-r.

If a needed primitive is missing upstream (e.g. Bessel functions, a constant),
prefer adding it upstream when it is general (dB / geometry), otherwise keep it
private in this crate.

### 2.2 Other dependencies
- `thiserror` 2 for `AntennaError` (validated constructors and parsing only).
- `rayon` (optional, feature `parallel`) for parallel batch evaluation.
- Dev: `approx`, `criterion` (benches).
- **No** `num-complex`: the array factor needs only `(re, im)` accumulation;
  use a tiny private `Complex` in `math.rs`. No special-function crate: write
  `bessel_j1` (and `j0` if needed) with documented accuracy (e.g. rational
  approximations, Abramowitz & Stegun 9.4) and test against tabulated values.
  Any new runtime dependency needs justification in an ADR.

### 2.3 Quality gates
- Follow the `rust-development`, `unit-testing-coverage`, `clean-code` skills
  on every code change; `secure-development` when touching dependencies,
  parsing of untrusted files (GRD) or CI; `application-documentation` for
  README and `docs/`; `requirements-analysis` when scoping a new model.
- `[lints.rust] unsafe_code = "forbid"`, `missing_docs = "warn"`;
  clippy `pedantic` (allowed with reasons in Cargo.toml: `module_name_repetitions`,
  index-to-f64 casts, `many_single_char_names`/`similar_names` for textbook notation).
- Coverage gate is 100 % *lines* across unit and integration tests. llvm-cov
  counts generic code once per test binary: keep non-hot helpers non-generic
  (`&dyn GainPattern` in `analysis`, slice-based parsers in `grd::ticra`) and
  make `should_panic` tests use a type that the positive tests also use.
- Release and bench profiles: `lto = "fat"`, `codegen-units = 1`.
- Before saying any work is done: `cargo fmt --all`, then `make check`
  (fmt, clippy `-D warnings`, tests with all features and default features,
  rustdoc `-D warnings`, 100 % line coverage via `cargo llvm-cov`,
  `cargo audit`, MSRV build, docs checks). Copy the Makefile, CI workflow
  (`.github/workflows/ci.yml`, Dependabot) and `scripts/docs/` from `rf-r` /
  `geometry-r` and adapt them.
- Documentation like the sibling crates: `README.md` (doctested via a
  `ReadmeDoctests` item in `lib.rs`), `docs/arc42/architecture.md`,
  `docs/c4/*.md` (Mermaid), `docs/technical-design.md`, `docs/adr/`
  (`scripts/docs/new_adr.sh "title"`), `CHANGELOG.md`. When a module or public
  type is added, removed or renamed, update the C4 component diagram and the
  technical-design component table.

### 2.4 Commands
```bash
cargo build
cargo test --all-features
cargo test <name>                                         # single test
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo doc --no-deps --open
cargo bench                                               # criterion (benches/)
make check                                                # everything CI runs
make msrv                                                 # cargo +1.88 check --lib
```

---

## 3. Conventions

### 3.1 Antenna frame and directions
- Right-handed antenna frame: **boresight = +Z**, reference axis +X, +Y = Z × X.
- θ = off-boresight angle in [0, π]; φ = azimuth in the X-Y plane from +X
  towards +Y, in [0, 2π). Direction cosines `u = sinθ cosφ`, `v = sinθ sinφ`,
  `w = cosθ`. Same as `geometry_r::ScanAngles`.
- Patterns are defined over the **whole sphere** (back hemisphere included).
  A model that is only defined in the front hemisphere must return its
  documented back-lobe level, never `NaN`.
- Principal cuts: φ = 0 is the "E/X cut" (X-Z plane), φ = π/2 the "H/Y cut"
  (Y-Z plane). Document per model which physical plane that is.
- Alternative grid conventions (az-over-el, el-over-az, θ/φ, u/v) are
  conversions on `Direction`, used by the GRD reader and by steering.

### 3.2 Gain
- All public gains are **absolute, in dBi**, including losses the model
  declares (efficiency, ohmic loss). A model exposes `peak_gain_dbi()`.
- "Relative gain" / "pattern" values are in dB relative to the peak (≤ 0) and
  are suffixed `_rel_db`.
- Losses are positive dB and subtracted (rf-r convention).
- Linear gain uses the `_lin` suffix; field amplitudes `_amp`.
- Floor: every pattern has a configurable minimum gain (`floor_dbi`) so nulls
  never return `-inf` unless the caller asks for that (`floor_dbi = -inf`).

### 3.3 Polarisation (scope boundary)
- Milestones M1–M5 compute **total gain** (co + cross power) unless the model
  states co-polar only. The data model for GRD and arrays keeps two complex
  field components so co/cross-polar gain (Ludwig-3, RHCP/LHCP) can be added
  later without breaking the API (§9).

### 3.4 API style
- Constructors validate and return `Result<_, AntennaError>`; evaluation
  methods take plain `f64`/`Direction`, are `#[inline]` + `#[must_use]`,
  never allocate, never panic on finite input.
- Non-finite directions propagate `NaN` (documented), they do not panic.
- Builders for models with many parameters (phased arrays, steering limits).
- No `unwrap()`/`expect()` outside tests. Batch functions document
  `# Panics` for slice-length mismatch.
- Every public item has `///` docs with the formula, its source
  (paper, ITU-R recommendation and version) and a runnable doc-test.

---

## 4. Functional requirements

IDs are stable. Reference them in tests and commits (e.g. `ANT-PAR-04`).

### 4.1 Core — `direction`, `pattern`, `error`, `math`
| ID | Requirement |
|----|-------------|
| ANT-CORE-01 | `Direction`: `Copy` unit vector in the antenna frame with constructors `from_theta_phi`, `from_uv` (front or back hemisphere), `from_vec3`, `from_az_el` / `from_el_az` (documented conventions), `from_tilts`; accessors for θ, φ, u, v, w. Lossless conversion to/from `geometry_r::ScanAngles`. |
| ANT-CORE-02 | `GainPattern` trait (object-safe, `Send + Sync`): `gain_dbi(&self, dir: Direction) -> f64`, `peak_gain_dbi(&self) -> f64`, and provided `gain_dbi_batch(&self, dirs: &[Direction], out: &mut [f64])`. Frequency is fixed at construction (no `frequency_hz` on the trait, ADR 0002); `PhasedArray::gain_dbi_at` covers wideband evaluation. Blanket impls for `&P`, `Box<P>`, `Arc<P>`. |
| ANT-CORE-03 | Static dispatch path: `AnyPattern` enum over all built-in models, implementing `GainPattern`, for callers who want to avoid `dyn`. Generic wrappers (`Steered<P>`, `PhasedArray<E>`) are generic over `P: GainPattern`. |
| ANT-CORE-04 | `AntennaError` (thiserror): `InvalidParameter { name, reason }`, `Parse { line, reason }`, `UnsupportedFormat`, `OutOfRange`, `Io` (only in `std::io` adapters). Tests assert specific variants. |
| ANT-CORE-05 | Pattern analysis utilities (any `GainPattern`): principal and arbitrary-φ cuts, gain map on a θ/φ or u/v grid (into caller buffers), numerical −3 dB / −10 dB beamwidth search, first-null and peak-sidelobe search on a cut. |
| ANT-CORE-06 | Numerical **directivity integration** over the sphere (`∫ G dΩ / 4π`), with a documented quadrature (Gauss–Legendre in θ, uniform in φ). Used to verify normalisation of every model and to normalise GRD data and arrays on demand. Isotropic → 0 dBi within 0.01 dB. |
| ANT-CORE-07 | `Isotropic { gain_dbi }` and `CosineQ { q, peak_gain_dbi }` (cos^q θ element, peak set by `2(q+1)` when not given) reference patterns. |

### 4.2 Parabolic / aperture antennas — `parabolic`
Peak gain uses `rf_r::antenna::parabolic_gain_dbi(D, f, η)`.

| ID | Requirement |
|----|-------------|
| ANT-PAR-01 | **Uniform circular aperture**: `G(θ) = G0 · [2 J1(x)/x]²`, `x = π D sinθ / λ`. Exact peak at θ = 0 (limit, no 0/0). |
| ANT-PAR-02 | **Tapered circular aperture**: parabolic-on-pedestal taper `(1 − (r/a)²)^n` with edge taper `C` (pedestal), n ∈ {0, 1, 2}; closed form via `Λ_n` functions. Aperture (taper) efficiency computed, combined with a user spill-over/ohmic efficiency. |
| ANT-PAR-03 | **Gaussian main beam**: `G = Gmax − 12 (θ/θ3dB)²` with a sidelobe floor; θ3dB from the user or from `k·λ/D` (k configurable, default 70°). |
| ANT-PAR-04 | **ITU-R reference patterns / envelopes** (one type each, parameters as in the recommendation, peak gain from D/λ or given), transcribed from the **latest edition in force** (checked on itu.int, 2026-11): **S.465-6** (Notes 4, 5), **S.580-6**, **S.1428-1**, **F.699-9**, **F.1245-3** (average + Annex 1 generalised), **S.672-4** recommends 1 (LN −20/−25), **S.1528-1** (recommends 1.2, 1.3 MEO/LEO, 1.4 Taylor), **BO.1213-1**, **BO.1295-0**, and **RR Appendices 30/30A** Edition 2024 (APS30: AP30 Annex 5 Figs. 7, 8, 9, 10, 11A/B, 13; AP30A Annex 3 Figs. A, B, C). Co/cross-polar components via `PolarComponent` where the text defines both. Below φmin, where S.465/S.580 define nothing, the RR App. 7 parabola floored at the envelope value is used (documented crate choice). Not modelled: S.672-4 recommends 2 (shaped beams, needs the coverage contour). Never transcribe from memory; test at and around each breakpoint and against the worked examples of the texts. |
| ANT-PAR-05 | **Elliptical beams**: `GaussianBeam::elliptical` and `ItuS672::elliptical` use an effective beamwidth per plane φ (ellipse with major/minor axes and orientation). No generic elliptical wrapper. |
| ANT-PAR-06 | Offset-fed / shaped reflectors are **not** modelled analytically: use a GRD file (§4.3). |

### 4.3 Custom gain from GRD files — `grd`
Target format: **TICRA GRASP `.grd` far-field grid** (ASCII), implemented from
TICRA's "Field data in rectangular grid" description (record layout, `IGRID`
and `ICOMP` definitions quoted in `grd/ticra.rs` and ADR 0003). Confirm other
"GRD" flavours with the owner before supporting them (§10).

| ID | Requirement |
|----|-------------|
| ANT-GRD-01 | Parse the TICRA ASCII grid: free-text header up to the `++++` line, then `KTYPE` (must be 1), `NSET ICOMP NCOMP IGRID`, per-set centres `IX IY`, then per set `XS YS XE YE`, `NX NY KLIMIT` (with `IS IN` row limits when `KLIMIT = 1`), then `NCOMP` complex components per point. Support multiple sets (beams / frequencies), selected by index. |
| ANT-GRD-02 | Grid types (`IGRID`): u/v (1), elevation-over-azimuth (4), elevation-and-azimuth (5), azimuth-over-elevation (6), θ/φ (7, `(X, Y) = (φ, θ)`). Each grid type maps a `Direction` to grid coordinates (x, y) through `direction` conversions; angles in degrees. |
| ANT-GRD-03 | Field components (`ICOMP`): linear θ/φ (1), RHCP/LHCP (2), linear co/cx Ludwig-3 (3) at minimum; total gain = `|E1|² + |E2|²`. Keep both components for later co/cross-polar support. Field normalisation follows GRASP (`10·log10|E|²` = dBi when power-normalised to 4π); document and allow a user offset or **re-normalisation by directivity integration** (ANT-CORE-06). |
| ANT-GRD-04 | Build a `GridPattern` from the parsed set: precompute gain (dB) on the grid, bilinear interpolation in grid coordinates, interpolation done in dB by default (configurable: dB or linear power). Optional bicubic later. Constant-time cell lookup (regular grid), no allocation at evaluation. |
| ANT-GRD-05 | Outside the grid: configurable policy `Floor(dBi)` (default: grid minimum), `Clamp` (nearest edge), or `NaN`. Periodicity in φ / azimuth is handled (no seam at 0/2π). |
| ANT-GRD-06 | Robust parsing of untrusted files: no panics on malformed input, explicit size limits (NX·NY·NSET), line numbers in `AntennaError::Parse`, Fortran-style numbers (`1.0D+01`, missing exponent letter). Fuzz-style tests on truncated/garbled files. |
| ANT-GRD-07 | Writer: export any `GainPattern` to a GRD grid (θ/φ or u/v) so round-trip tests and external tools can consume crate patterns. |
| ANT-GRD-08 | Simple tabular input as a second source: CSV of (θ_deg, φ_deg, gain_dBi) on a regular grid, and 1-D cuts (θ, gain) assumed rotationally symmetric. Same `GridPattern` result. |

### 4.4 Steerable antennas — `steer`
Wraps any pattern `P` and a pointing; gain is evaluated by rotating the
requested direction from the **mount frame** into the pattern frame.

| ID | Requirement |
|----|-------------|
| ANT-STR-01 | `Steered<P>`: pattern + mount-to-antenna rotation (`geometry_r::Mat3`, precomputed). Pointing setters: by `Direction` (boresight target), az/el, X/Y tilt, θ/φ, plus a rotation about boresight (polarisation / pattern roll). Changing pointing is O(1) and allocation-free. |
| ANT-STR-02 | Mechanical mount types: az/el (keyhole at zenith documented), X/Y (elevation-over-cross-elevation), hexapod/free (direct rotation). Each mount converts a target direction into axis angles and back. |
| ANT-STR-03 | Steering limits: per-axis angle ranges and a maximum scan cone. Pointing outside limits returns `AntennaError::OutOfRange` from the validating setter; an unchecked fast setter exists and clamps. |
| ANT-STR-04 | **Electronic scan loss** option for electronically steered antennas modelled as a single pattern: gain reduced by `−10·n·log10(cos θ_scan)` (n configurable, default 1, i.e. projected-aperture loss) and beam broadening `θ3dB / cos θ_scan` in the scan plane. Max scan angle enforced. |
| ANT-STR-05 | Mispointing: a fixed pointing error (bias) and helpers to evaluate gain loss for a given error, for tracking-error studies. |

### 4.5 Phased arrays — `array`
`G(û) = G_e(û) · |AF(û)|² / Σ|w_n|²`, with `AF(û) = Σ w_n · exp(j k r_n · û)`,
`k = 2π/λ`. This normalisation assumes the element pattern `G_e` is the
embedded (in-array) element gain and holds for well-sampled lattices without
mutual coupling; document it, and provide exact directivity normalisation via
ANT-CORE-06 as an option.

| ID | Requirement |
|----|-------------|
| ANT-PAA-01 | `PhasedArray<E: GainPattern>`: element positions (`Vec3`, metres, antenna frame), one shared element pattern (optionally per-element orientation later), complex weights. Validated by a builder. |
| ANT-PAA-02 | Lattice generators: linear, rectangular, triangular/hexagonal, concentric rings, with spacing in metres or in λ at a design frequency; circular or arbitrary-polygon aperture boundary (element kept if inside). Custom positions accepted. |
| ANT-PAA-03 | Amplitude tapers: uniform, cosine-on-pedestal, **Taylor** (n̄, SLL), **Dolph–Chebyshev** (linear; separable rectangular via `separable`); user weights. Report taper efficiency. Circular Taylor moved to the roadmap. |
| ANT-PAA-04 | Steering: phase steering (phases computed at a design frequency → **beam squint** when evaluated at another frequency) and true-time-delay steering (no squint). Steering by `Direction`. Optional **phase quantisation** (N bits) and amplitude quantisation. |
| ANT-PAA-05 | Evaluation at an arbitrary frequency `f` (wideband behaviour): `gain_dbi_at(dir, freq_hz)`; `gain_dbi(dir)` uses the design frequency. |
| ANT-PAA-06 | Impairments: element failures (mask), random amplitude/phase errors from a seeded RNG supplied by the caller (no global RNG, no new dependency: accept pre-generated error vectors). |
| ANT-PAA-07 | Analysis helpers: grating-lobe check for a lattice/spacing/scan range (visible-space test of `û₀ + λ·lattice reciprocal`), scan loss vs. angle, broadside directivity estimate `4πA/λ²`. |
| ANT-PAA-08 | Multiple simultaneous beams: `MultiBeamArray` sharing positions, one weight set per beam; evaluate one beam or all beams per direction. |
| ANT-PAA-09 | Performance: precompute `k·r_n` once per frequency; inner loop is a branch-free `sin_cos` accumulation over contiguous `x, y, z, w_re, w_im` arrays (SoA layout); batch over directions; `parallel` feature parallelises over directions. Benched in `benches/antenna_bench.rs` (1024 elements). The separable fast path (`AF = AF_x · AF_y`) moved to the roadmap. |

### 4.6 Integration with geometry — `frames`
| ID | Requirement |
|----|-------------|
| ANT-FRM-01 | Gain towards an ECEF target from a `geometry_r::SatelliteFrame` (body frame = antenna frame, or with a fixed antenna-to-body mounting rotation). |
| ANT-FRM-02 | Gain towards a satellite from a ground station `geometry_r::coords::LocalFrame` (ENU) with a steerable antenna pointed at another target (off-axis gain for interference, e.g. adjacent satellite). |
| ANT-FRM-03 | Batch versions over many targets, reusing precomputed frames. |

---

## 5. Non-functional requirements
| ID | Requirement |
|----|-------------|
| ANT-NFR-01 | Evaluation of analytic patterns (parabolic, ITU-R, grid lookup): no allocation, O(1); bench-backed. |
| ANT-NFR-02 | Thread safety: all models `Send + Sync`, immutable after build; steering mutation via `&mut self` only. |
| ANT-NFR-03 | Determinism: identical results for identical inputs across runs and with/without `parallel` (bitwise for per-direction work). |
| ANT-NFR-04 | Accuracy targets: analytic patterns within 1e-9 dB of the reference formula; Bessel J1 within 1e-8 absolute; GRD interpolation exact at grid nodes. |
| ANT-NFR-05 | Every optimisation is backed by a criterion benchmark in `benches/`. |

---

## 6. Architecture

```
src/
  lib.rs           crate docs, module list, prelude, ReadmeDoctests
  error.rs         AntennaError, Result
  math.rs          private Complex, bessel_j0/j1, quadrature, interpolation
  direction.rs     Direction, grid conventions (θφ, uv, az/el, el/az, tilts), ScanAngles interop
  pattern.rs       GainPattern trait, AnyPattern, Isotropic, CosineQ
  analysis.rs      cuts, gain maps, beamwidth / null / sidelobe search, directivity integration
  parabolic/
    mod.rs         UniformAperture, TaperedAperture, GaussianBeam, elliptical mapping
    itu/           ITU-R patterns: mod.rs (shared envelopes, PolarComponent), earth.rs (S.465,
                   S.580, S.1428), fixed.rs (F.699, F.1245), satellite.rs (S.672, S.1528),
                   bss.rs (BO.1213, BO.1295), rr_ap30.rs (RR Appendices 30/30A)
  grd/
    mod.rs         GridPattern, interpolation, out-of-grid policy
    ticra.rs       TICRA GRASP .grd reader / writer
    table.rs       CSV / 1-D cut reader
  steer.rs         Steered<P>, mounts, limits, scan loss, mispointing
  array/
    mod.rs         PhasedArray<E>, MultiBeamArray, SteeringMode, evaluation kernel
    lattice.rs     lattice generators and aperture masks
    taper.rs       uniform, cosine-on-pedestal, Taylor, Chebyshev, separable, power efficiency
    grating.rs     grating-lobe prediction
    excitation.rs  hardware limits: phase-only, attenuator range / bits, phase quantisation
    synthesis.rs   pattern synthesis (alternating projections), closed-form nulls
  frames.rs        SatelliteFrame / LocalFrame integration
  batch.rs         par::* (feature `parallel`)
benches/           criterion benches per family
tests/             integration scenarios (§7)
examples/          dish_pattern.rs, grd_lookup.rs, steered_dish.rs, phased_array_scan.rs
docs/              arc42/, c4/, technical-design.md, adr/
```

ADRs (in `docs/adr/`):
1. Dependency on `rf-r` and `geometry_r` (git deps, MSRV 1.88 inherited).
2. `GainPattern` trait + `AnyPattern` enum (dyn and static dispatch).
3. Antenna-frame and direction conventions (aligned with `ScanAngles` and TICRA).
4. Array gain normalisation (`|AF|²/Σ|w|²` × embedded element gain).
5. GRD interpolation in dB and out-of-grid policy.
6. ITU-R patterns from the latest editions.
7. Array pattern synthesis by alternating projections.

---

## 7. Testing requirements
- Unit tests at the bottom of each module (`#[cfg(test)] mod tests`), shared
  `assert_close(actual, expected, tol)` helper; integration tests in `tests/`.
  Cover positive, negative and boundary cases for every public function.
- Mandatory edge cases: θ = 0 (0/0 limits), θ = π (back direction), φ wrap at
  0/2π, u² + v² = 1 (horizon), non-finite inputs, zero/negative diameters,
  frequencies and efficiencies, empty arrays, single-element arrays, slice
  length mismatch, grid edges and corners, malformed GRD files.
- Reference values and identities:
  - Uniform circular aperture: first null at `sinθ = 1.2197 λ/D`, first
    sidelobe −17.57 dB, −3 dB full beamwidth ≈ `1.029 λ/D` rad; peak equals
    `rf_r::antenna::parabolic_gain_dbi`.
  - Directivity integration: isotropic 0 dBi, cos^q element `10·log10(2(q+1))`,
    every analytic model integrates to its declared peak/efficiency within
    tolerance.
  - ITU-R envelopes: values at and around each breakpoint, hand-computed from
    the cited recommendation text; continuity where the recommendation says so.
  - Arrays: N-element uniform linear array at λ/2 → broadside directivity
    `10·log10 N` (isotropic elements); Chebyshev design → all sidelobes at the
    design level; phase steering at f ≠ f0 → squint `sinθ = (f0/f) sinθ0`;
    grating lobe appears at spacing `d > λ/(1+|sinθ0|)`; separable fast path
    equals the general kernel.
  - Steering: gain at the pointed direction equals the peak; rotation round
    trips (target → axis angles → target); scan loss `cos^n`.
  - GRD: synthetic file written by ANT-GRD-07 round-trips; values exact at
    nodes; bilinear midpoint check; every IGRID/ICOMP combination supported.
- Property-style checks: `Direction` conversions round-trip through every
  convention and through `ScanAngles`; batch == scalar; `parallel` == serial.

---

## 8. Milestones
1. **M0 Scaffolding**: done.
2. **M1 Parabolic**: done.
3. **M2 ITU-R envelopes**: done for every pattern listed in ANT-PAR-04 (latest
   editions, verified against the official texts).
4. **M3 GRD**: done (ICOMP 1–4, IGRID 1/4/5/6/7, CSV and cuts, writer).
5. **M4 Steerable**: done, with frames integration.
6. **M5 Phased arrays**: done except circular Taylor and the separable fast path.
7. **M6 Performance**: criterion benches for every family and the `parallel`
   feature done; allocation audit of evaluation paths done (no allocation).

---

## 9. Roadmap (after M6)
- Co-polar / cross-polar gain (Ludwig-3, RHCP/LHCP), axial ratio, polarisation
  mismatch loss — using the two field components kept since M3/M5.
- Per-element orientation and per-element patterns (GRD per element) in arrays.
- Bicubic / spline GRD interpolation; TICRA `.cut` files; GRD `ICOMP` 5–9.
- Circular Taylor taper; separable fast path for rectangular arrays.
- More ITU-R patterns on demand (e.g. S.672-4 recommends 2 shaped beams, F.1336,
  S.1855, BO.1443), from verified text only; re-check editions on itu.int.
- `serde` support for model configurations (feature).
- Sub-arrays, overlapped sub-arrays, digital beamforming weight solvers
  (MVDR). Done: per-element excitation, hardware limits, mask synthesis by
  alternating projections, closed-form null steering (ADR 0007).
- Frequency-dependent analytic patterns (gain vs. frequency over a band).

---

## 10. Open questions for the owner
- "GRD" format: TICRA GRASP grid (implemented) or also another tool's export
  (e.g. a measurement-range format, Surfer `.grd`)? The GRASP 4π W field
  normalisation was confirmed by the owner (2026-11-07).
- Crate/package name: `antenna_r` (assumed, like `geometry_r`) or `antenna-r` (like `rf-r`)?
- Should `geometry_r` stay a mandatory dependency (it pulls `i_overlay` and
  `serde_json`), or be behind a default feature with a minimal local `Vec3`?
  Assumed mandatory for now (ADR 1).
