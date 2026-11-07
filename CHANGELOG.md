# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses semantic versioning.

## [Unreleased]

### Added
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
