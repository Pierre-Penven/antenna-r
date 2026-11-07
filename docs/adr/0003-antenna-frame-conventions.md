# 0003. Antenna-frame and grid conventions

Date: 2026-11-06
Status: Accepted

## Context

Patterns, pattern files, positioners and satellite frames all describe directions. Mismatched
conventions silently mirror or rotate patterns.

## Decision

- Antenna frame: right-handed, boresight `+Z`, reference `+X`. `θ` from `+Z`, `φ` from `+X`
  towards `+Y`; `u = sinθ cosφ`, `v = sinθ sinφ`. Identical to `geometry_r::ScanAngles`.
- GRD grid coordinates follow the TICRA GRASP definitions verbatim (from TICRA's "Field data in
  rectangular grid" description): `IGRID 1` (u, v); `4` elevation over azimuth
  `r̂ = (−sin Az cos El, sin El, cos Az cos El)`; `5` `Az = −θ cos φ, El = θ sin φ`;
  `6` azimuth over elevation `r̂ = (−sin Az, cos Az sin El, cos Az cos El)`; `7` `(X, Y) = (φ, θ)`.
  Angles in the files are in degrees.
- Steering mounts: az/el uses compass azimuth (clockwise from `+Y` towards `+X`) so that an ENU
  mount frame gives the usual azimuth from north, matching `geometry_r::LocalFrame`.

## Consequences

- `Direction` converts losslessly to and from `ScanAngles`, so satellite scan angles can be fed
  to any pattern directly.
- GRD files exported by GRASP are read without axis flips; θ ∈ [−180°, 180°] grids are handled
  through the `(−θ, φ + 180°)` identity.
