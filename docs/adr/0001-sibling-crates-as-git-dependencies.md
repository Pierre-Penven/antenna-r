# 0001. Sibling crates as pinned git dependencies

Date: 2026-11-06
Status: Accepted

## Context

Antenna patterns need dB math and aperture gain (in `rf-r`) and vectors, rotations, scan angles and
satellite/ground frames (in `geometry_r`). Neither crate is on crates.io yet.

## Options considered

1. Copy the needed code into this crate: no dependency, but duplicated, diverging code.
2. Path dependencies on local checkouts: does not build for anyone else or in CI.
3. Git dependencies pinned to a revision.

## Decision

Option 3: `rf-r` and `geometry_r` as git dependencies with a pinned `rev`. `geometry_r` is a
mandatory dependency (not behind a feature) so that `Direction` can reuse `Vec3`, `Mat3` and
`ScanAngles` directly.

## Consequences

- No duplicated physics or geometry code; conventions match the sibling crates by construction.
- The MSRV is inherited from `geometry_r`: Rust 1.88.
- `geometry_r` brings `i_overlay` and `serde_json` into the build, unused here.
- Switch to crates.io versions when the siblings are published; bump the pinned revisions explicitly.
