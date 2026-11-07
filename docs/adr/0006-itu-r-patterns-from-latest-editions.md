# 0006. ITU-R patterns from the latest editions

Date: 2026-11-07
Status: Accepted

## Context

The owner asked for every ITU-R pattern on the roadmap, in its latest version. Several texts were
revised recently (S.1528-1 in 08/2026, F.699-9 in 02/2026), the PDFs render some formulas badly
(lost signs, square roots and Greek letters), and some recommendations split the pattern into
co-polar and cross-polar components.

## Options considered

1. Transcribe from secondary sources or memory: fast, but unverifiable.
2. Transcribe from the official PDFs of the editions in force, check garbled formulas on the
   rendered pages and against the worked examples and continuity of each text.

## Decision

Option 2. Implemented editions: S.465-6, S.580-6, S.1428-1, F.699-9, F.1245-3, S.672-4,
S.1528-1, BO.1213-1, BO.1295-0, and the Radio Regulations (Edition 2024) Appendix 30 Annex 5 and
Appendix 30A Annex 3 patterns. Each type cites its edition. Patterns defined per polarisation take
a `PolarComponent`; `peak_gain_dbi` returns the maximum of the selected component. Where a text is
ambiguous the crate documents its choice:

- S.1528-1 `recommends 1.4`: `σ` uses `μl` (matches the Annex 2 example), not the `μ1` printed.
- S.1528-1 `recommends 1.3`: the `Gm` form of the side lobes is used; the `D/λ` form coincides when
  `Gm ≈ 20 log(D/λ) + 7.97`.
- S.1428-1 leaves 33.1° undefined for `20 ≤ D/λ ≤ 25`: −9 dBi is used.

Checking the texts again exposed a missing square root in S.672-4 Table 1
(`a = 2.58 √(1 − log z)`), fixed in the same change; circular beams were not affected.

## Consequences

- Each pattern is tested at its breakpoints and against the examples of its text (S.1528 Annex 1
  and 2, BO.1213 60 cm and 45 cm antennas, F.1245 Note 1).
- New editions must be re-checked on itu.int before a release; the editions are listed in the
  module documentation of `parabolic::itu`.
- S.672-4 `recommends 2` (shaped beams) needs a coverage contour and stays out of scope.
