# 0005. GRD interpolation and out-of-grid policy

Date: 2026-11-06
Status: Accepted

## Context

Gridded patterns must be evaluated between nodes and outside their coverage, quickly and safely.

## Decision

- Bilinear interpolation in grid coordinates, on gain in dB by default (`BilinearLinear` on linear
  power is available). Values are exact at nodes. Zero fields are stored as −300 dBi so dB
  interpolation stays finite; missing points (`KLIMIT = 1`) are `NaN` and a missing corner with a
  non-zero weight makes the point "outside".
- Azimuth-like axes wrap every 360°, and a full-turn axis closes its last cell onto the first.
- Outside policy: `GridMinimum` (default), `Floor(dBi)`, `Clamp` (nearest edge) or `Nan`.
- Field normalisation: `10 log10(|F1|² + |F2|²)` is taken as dBi (GRASP convention for 4π W);
  `with_offset_db` and `normalised_to_directivity` cover other normalisations.

## Consequences

- O(1) lookup with no allocation.
- Bilinear dB interpolation slightly misrepresents very deep nulls; bicubic interpolation is on
  the roadmap.
