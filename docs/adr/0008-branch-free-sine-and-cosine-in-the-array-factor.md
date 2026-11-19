# 0008. Branch-free sine and cosine in the array factor

Date: 2026-11-19
Status: Accepted

## Context

An exact array gain sums one term per element, `a_n exp(jφ_n)`, and the library sine and cosine
(`f64::sin_cos`) take most of that time: they are calls with branches, so the loop runs one element
at a time. Users that need many gains (an active-array payload model asks for tens of thousands per
evaluation, each over hundreds of elements) are bound by this sum.

## Options considered

1. Keep `f64::sin_cos`.
2. A SIMD crate or explicit intrinsics: a new dependency or `unsafe`, both excluded by `CLAUDE.md`.
3. Our own sine and cosine, written without branches or calls, and the sum split into independent
   partial sums over fixed-size chunks so that the compiler vectorises it on the baseline target.
4. The shift-invariant kernel (`array::kernel`): table lookups, but approximate and built per
   frequency; it stays available for uses that accept that.

## Decision

Option 3. `math::branch_free_sin_cos` uses the three-part Cody–Waite reduction by π/2 of fdlibm
(`__ieee754_rem_pio2`, exact products for `round(2x/π) < 2^20`), the fdlibm minimax polynomials
`__kernel_sin` / `__kernel_cos` on [−π/4, π/4] and the quadrant as bit masks. It is within `1e-15`
of `f64::sin_cos` for `|x| ≤ 10⁶` rad (tested over quadrant boundaries and 200 000 magnitudes from
1e-9 to 10⁶). `PhasedArray::array_factor` sums in four partial sums over `as_chunks::<4>()`; if any
element phase exceeds 10⁶ rad (elements kilometres apart) it redoes the sum in order with
`f64::sin_cos`.

## Consequences

- About 3.9× faster on `phased_array_1024_elements_1k_dirs` (12.7 → 3.2 ms) and 3.3× on
  `array_400_exact_gain` (4.14 → 1.26 µs), x86-64 baseline target, no new dependency, no `unsafe`.
- Results change by rounding only: each term within `1e-15`, the sum in a different order (gains
  within about 1e-12 dB). They stay deterministic and identical with and without `parallel`.
- The sign of a zero sine is not kept (it cannot change a sum).
