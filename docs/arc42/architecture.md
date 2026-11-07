# antenna-r architecture (arc42)

## 1. Introduction and goals

`antenna_r` computes far-field antenna gain in any direction for parabolic and aperture antennas,
gridded (GRD) patterns, steerable antennas and phased arrays, behind one trait (`GainPattern`).
It is a building block for link budgets, interference studies and coverage simulations.

Quality goals, in priority order:

| Priority | Goal | Scenario |
|---|---|---|
| 1 | Correctness | Every formula is traceable to its source (ITU-R text, TICRA format description, textbook) and checked against reference values in tests. |
| 2 | Performance | Gain evaluation does not allocate and does not return `Result`; batch and parallel variants exist. |
| 3 | Composability | Any pattern can be steered or used as an array element; callers do not depend on the antenna type. |
| 4 | Robustness | Untrusted GRD/CSV input never panics and has bounded memory use. |

## 2. Constraints

- Rust edition 2024, MSRV 1.88 (inherited from `geometry_r`), `unsafe` forbidden.
- Library only: no I/O beyond parsing caller-supplied text, no async, no global state.
- Sibling crates `rf-r` and `geometry_r` are git dependencies pinned to a revision until they are
  published on crates.io ([ADR 0001](../adr/0001-sibling-crates-as-git-dependencies.md)).
- MIT licence; publishing to crates.io needs the owner's explicit go-ahead.

## 3. Context and scope

See the [C4 context diagram](../c4/1-context.md). The crate is linked into the caller's program.
It depends on `rf-r` (dB conversions, aperture gain, wavelength, speed of light) and `geometry_r`
(`Vec3`, `Mat3`, `ScanAngles`, `SatelliteFrame`, `LocalFrame`). Pattern files come from external
tools (TICRA GRASP, measurement ranges) as text passed in by the caller.

## 4. Solution strategy

- One object-safe trait, `GainPattern`, with absolute gain in dBi
  ([ADR 0002](../adr/0002-gain-pattern-trait-and-any-pattern.md)); `AnyPattern` for static dispatch.
- Construct once, evaluate many times: constructors validate and precompute (`Result`), evaluation
  is plain `f64` math, `#[inline]`, allocation-free.
- One antenna frame convention shared with `geometry_r::ScanAngles` and TICRA GRASP
  ([ADR 0003](../adr/0003-antenna-frame-conventions.md)).
- Formulas transcribed from primary sources only (latest ITU-R editions, Radio Regulations 2024),
  checked against the worked examples of each text ([ADR 0006](../adr/0006-itu-r-patterns-from-latest-editions.md)).

## 5. Building block view

See the [C4 container](../c4/2-containers.md) and [component](../c4/3-components.md) diagrams.

| Module | Responsibility |
|---|---|
| `direction` | `Direction` unit vector, θ/φ, u/v, tilts and TICRA az/el conventions |
| `pattern` | `GainPattern`, `AnyPattern`, `Isotropic`, `CosineQ` |
| `parabolic` | `UniformAperture`, `TaperedAperture`, `GaussianBeam` |
| `parabolic::itu` | ITU-R S.465, S.580, S.1428, F.699, F.1245, S.672, S.1528, BO.1213, BO.1295 and RR Appendix 30/30A patterns |
| `grd` | `GridPattern`, `GridType`, `GridAxis`; `ticra` reader/writer; `table` CSV and cuts |
| `steer` | `Steered<P>`, `Mount`, limits, `ScanLoss`, pointing error |
| `array` | `PhasedArray<E>`, `MultiBeamArray`; `lattice`, `taper`, `grating` |
| `frames` | gain towards ECEF targets from satellite or ground frames |
| `analysis` | cuts, gain maps, beamwidths, nulls, side lobes, power integration |
| `math` | `Complex`, Bessel `J_n`, `Λ_n`, Gauss–Legendre, Chebyshev |
| `batch` | rayon parallel batches (feature `parallel`) |
| `error` | `AntennaError` |

## 6. Runtime view

**Gain of a steered dish towards a satellite** (`frames::gain_from_station`):
1. `LocalFrame::to_enu` turns the ECEF target into an ENU vector → `Direction` (mount frame).
2. `Steered::gain_dbi` rotates it into the antenna frame with the precomputed matrix (`Rᵀ d`),
   applies scan-plane broadening if configured, and calls the inner pattern.
3. The inner pattern evaluates its closed form (e.g. `Λ_n` for a tapered aperture) and applies its floor.

**Reading a GRD file**: `GrdFile::parse` validates the header, counts and limits record by record,
then `GrdFile::pattern(i)` converts the two field components to gain in dB and builds a
`GridPattern`, which maps each `Direction` to grid coordinates and interpolates bilinearly.

**Phased array evaluation**: `PhasedArray::gain_dbi_at` sums `a_n exp(j(k r_n·û + ψ_n))` over
structure-of-arrays element data and adds the element gain
([ADR 0004](../adr/0004-array-gain-normalisation.md)).

## 7. Deployment view

None: the crate is compiled into the caller's binary. CI runs on GitHub Actions
(`.github/workflows/ci.yml`).

## 8. Crosscutting concepts

- **Units**: SI inside, unit suffix on every name (`_dbi`, `_db`, `_rad`, `_deg`, `_hz`, `_m`).
- **Errors**: one `AntennaError` enum (thiserror); only constructors, setters and parsers fail.
- **Non-finite input**: directions with `NaN` give `NaN` gains, never a panic.
- **Floors**: patterns with nulls return at least `floor_dbi` (default −100 dBi).
- **Input validation**: GRD parsing caps points (`MAX_GRID_POINTS`) and sets (`MAX_SETS`),
  checks every count and returns line numbers.
- **Testing**: unit tests per module, integration scenarios in `tests/`, README doctests,
  100 % line coverage enforced in CI.

## 9. Architecture decisions

- [ADR 0001](../adr/0001-sibling-crates-as-git-dependencies.md): sibling crates as pinned git dependencies.
- [ADR 0002](../adr/0002-gain-pattern-trait-and-any-pattern.md): `GainPattern` trait and `AnyPattern`.
- [ADR 0003](../adr/0003-antenna-frame-conventions.md): antenna-frame and grid conventions.
- [ADR 0004](../adr/0004-array-gain-normalisation.md): phased-array gain normalisation.
- [ADR 0005](../adr/0005-grd-interpolation-and-out-of-grid-policy.md): GRD interpolation and out-of-grid policy.
- [ADR 0006](../adr/0006-itu-r-patterns-from-latest-editions.md): ITU-R patterns from the latest editions.

## 10. Quality requirements

| Scenario | Target |
|---|---|
| Analytic pattern accuracy | within 1e-9 dB of the reference formula |
| Bessel functions | 1e-12 absolute |
| GRD interpolation | exact at grid nodes |
| Uniform aperture | first null 1.2197 λ/D, side lobe −17.57 dB, HPBW 1.029 λ/D |
| Array | λ/2 isotropic ULA directivity = N; Chebyshev side lobes equal to the design |
| Robustness | malformed GRD/CSV input returns `AntennaError`, never panics |

## 11. Risks and technical debt

- GRASP field normalisation (`|E|²` = gain for 4π W) is used (confirmed by the owner); files
  normalised differently need an offset or `normalised_to_directivity`.
- ITU-R texts are revised (S.1528-1 and F.699-9 date from 2026): editions must be re-checked
  before each release.
- Array normalisation ignores mutual coupling; use an embedded element pattern from a GRD file.
- No separable fast path for rectangular arrays yet; evaluation is O(N) per direction.

## 12. Glossary

| Term | Meaning |
|---|---|
| dBi | gain relative to an isotropic radiator |
| Boresight | antenna `+Z` axis, direction of maximum gain for a symmetric beam |
| u, v | direction cosines along antenna X and Y |
| HPBW | half-power (−3 dB) beamwidth |
| Taper | amplitude distribution over an aperture or array |
| Squint | beam pointing shift with frequency under phase steering |
| GRD | TICRA GRASP far-field grid file |
| TTD | true time delay steering |
