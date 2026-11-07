# 0002. `GainPattern` trait and `AnyPattern` enum

Date: 2026-11-06
Status: Accepted

## Context

Callers must evaluate very different antennas through one interface, compose them (steer any
pattern, use any pattern as an array element) and evaluate millions of directions.

## Options considered

1. A closed enum of all models: fast static dispatch, but no user models and no generic wrappers.
2. An object-safe trait only: open and composable, but `dyn` calls cannot be inlined.
3. A trait plus an enum of the built-in leaf models.

## Decision

Option 3. `GainPattern` (object-safe, `Send + Sync`) has `gain_dbi(Direction)`,
`peak_gain_dbi()` and a provided `gain_dbi_batch`. Wrappers are generic (`Steered<P>`,
`PhasedArray<E>`) so they inline the inner pattern. `AnyPattern` covers the built-in leaf models.
Blanket implementations for `&P`, `Box<P>` and `Arc<P>` allow trait objects anywhere.
Frequency is not part of the trait: models fix it at construction, and `PhasedArray` adds
`gain_dbi_at(dir, freq_hz)` for wideband evaluation. Analysis helpers take `&dyn GainPattern`:
they are not hot paths and a single compiled copy keeps them simple.

## Consequences

- User-defined patterns work with every wrapper and analysis function.
- Adding a built-in model means adding an `AnyPattern` variant.
