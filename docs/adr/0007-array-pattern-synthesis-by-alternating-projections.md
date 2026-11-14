# 0007. Array pattern synthesis by alternating projections

Date: 2026-11-14
Status: Accepted

## Context

Contoured (shaped) beams and side-lobe / not-to-exceed regions are specified as lower and upper
gain bounds over regions of directions, for arrays of any lattice (custom positions included)
up to a few thousand elements, under beamformer limits (phase-only, attenuator range and
resolution, phase-shifter resolution, failed elements). Lower bounds on `|F|` make the problem
non-convex.

## Options considered

1. Closed-form tapers (Taylor, Chebyshev) and Woodward–Lawson: fast, but only for regular
   lattices and simple shapes.
2. Convex optimisation (SOCP): globally optimal for upper bounds with a fixed main-beam phase,
   but lower (coverage) bounds are non-convex; needs a solver dependency.
3. Alternating projections (Bucci et al., IEE Proc. H 1990): project the field onto the bounds,
   project the excitation onto the hardware set, repeat. Any lattice and any bounds; local
   optimum; hardware limits are just another projection.
4. Global stochastic search (genetic, particle swarm): no guarantees, much slower.

## Decision

Option 3, matrix-free (fields and adjoints computed on the fly, rows in parallel), with a
least-squares back-projection of four conjugate-gradient steps (CGLS) per iteration — a single
gradient step converged far too slowly once element patterns make the problem ill-conditioned —
and the best iterate kept. Bounds can be relative to a coverage level, fixed or maximised:
shaping at the level fitted to the coverage samples (aimed at their power mean, the flat-top
level by power conservation, until the weakest sample reaches it), then fixed-level runs judged
at the level they reach: at and up to 8 dB below the shaped level when it misses the bounds,
then upward (1 dB steps, bisection to 0.1 dB) while they are met. Closed-form nulls (linearly
constrained minimum change) are a separate function. The result depends on the start: a pencil
beam far narrower than the coverage converges to a poor local optimum (a few dB below the flat
top), so `coverage_start` gives a matched start (least-squares fit of a flat, in-phase field over
the coverage, Woodward–Lawson style) for callers without a better excitation.

## Consequences

- Memory O(samples + elements); cost O(samples × elements) per iteration.
- Results depend on the start (warm start from the array's current excitation) and are local
  optima; reproducible (no randomness).
- No new dependency.
