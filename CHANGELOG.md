# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses semantic versioning.

## [Unreleased]

### Added
- `array::kernel`: `PhasedArray::kernel` and `ArrayKernel` — the shift-invariant gain of a planar
  array without phase quantisation (natural or shaped beams, phase or true-time-delay steering,
  squint off the design frequency), tabulated once as two 1D tables (rectangular grid with
  row × column weights, the separable fast path) or one 2D table; a gain costs a few lookups
  (≈ 12 ns against ≈ 6 µs for 400 elements).
- Array excitation: `PhasedArray::with_excitation` (per-element amplitude and phase offset added
  to the steering phase), accessors (`excitation`, `positions`, `failed`, `phase_bits`,
  `steering_phases_rad`, `radiated_amplitudes`), `power_efficiency` (per-element amplifiers).
- `array::excitation`: `ExcitationLimits` (phase-only, attenuator range and resolution),
  `quantise_phase`.
- `array::synthesis`: `synthesize` (lower / upper gain bounds at sampled directions, coverage level
  fixed or maximised, hardware limits, warm start) and `place_nulls` (closed-form nulls).
- `array::synthesis::coverage_start`: matched start for a coverage far wider than the beam
  (least-squares fit of a flat, in-phase field over the lower-bounded samples).
- `analysis::find_peak` and `analysis::refine_peak` (true pattern maximum).
- `math::Complex`: `ZERO`, `conj`, `scale`, `abs`, subtraction.

### Changed
- Synthesis with a maximised coverage level aims at the power mean of the coverage gains (was the
  median) and keeps shaping until the weakest sample reaches it: a coverage-only goal far wider
  than the start beam is now filled instead of being returned unchanged as "feasible".
- `PhasedArray::with_phase_bits` quantises the total phase (steering + excitation offset); with
  zero offsets the result is unchanged.

- `GainPattern` trait, `AnyPattern`, `Direction` with TICRA and `ScanAngles` conventions.
- Parabolic models: uniform and tapered circular apertures, circular and elliptical Gaussian beams.
- ITU-R S.465-6 (with Notes 4 and 5), S.580-6 and S.672-4 (recommends 1) patterns.
- ITU-R S.1428-1, F.699-9, F.1245-3 (average and Annex 1 generalised), S.1528-1 (recommends 1.2,
  1.3 and 1.4), BO.1213-1 and BO.1295-0 patterns, and the Radio Regulations Appendix 30/30A
  (APS30) earth-station and satellite patterns, with `PolarComponent` for co/cross-polar gain.
- TICRA GRASP `.grd` reader and writer, CSV θ/φ tables and 1-D cuts, `GridPattern` interpolation.
- Steerable antennas: free, az/el and X-over-Y mounts, limits, roll, pointing error, scan loss.
- Phased arrays: lattices, uniform/cosine/Taylor/Chebyshev tapers, phase and TTD steering,
  phase quantisation, failures and errors, multi-beam, grating-lobe prediction.
- Frames integration with `geometry_r`, pattern analysis helpers, optional `parallel` feature.

### Fixed
- ITU-R S.672-4 elliptical beams: Table 1 `a` now includes its square root
  (`2.58 √(1 − log z)`); circular beams were unaffected.
