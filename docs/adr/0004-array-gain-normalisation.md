# 0004. Phased-array gain normalisation

Date: 2026-11-06
Status: Accepted

## Context

Array gain needs a normalisation of `|AF|²`. Exact directivity requires integrating over the
sphere for every steering, which is too slow for routine use.

## Options considered

1. Numerical directivity per steering state.
2. `G = G_e · |AF|² / Σ|a_n|²` with `G_e` the embedded element gain.

## Decision

Option 2 by default. It equals the exact gain for well-sampled lattices (spacing ≤ λ/2 or no
grating lobes) when `G_e` is the embedded element gain, and it is what the λ/2 ULA test checks
(directivity `N`). Failed elements are excluded from `Σ|a_n|²`: gain is relative to the power fed
to radiating elements. `analysis::directivity_dbi` remains available for exact checks.

## Consequences

- O(N) per direction, no integration.
- Mutual coupling is not modelled; supply a measured or simulated embedded element pattern
  (e.g. a `GridPattern`) to account for it.
