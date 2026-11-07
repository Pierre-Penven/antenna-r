# C4 level 1: System context

`antenna_r` is a library linked into the developer's program. It has no runtime connections to
other systems; pattern files reach it as text read by the caller.

```mermaid
C4Context
  title System Context: antenna_r
  Person(developer, "Application developer", "Builds link-budget, interference or coverage tools in Rust")
  System(antenna, "antenna_r", "Rust library: antenna gain patterns behind one trait")
  System_Ext(rfr, "rf-r", "dB math, aperture gain, wavelength")
  System_Ext(geo, "geometry_r", "Vectors, rotations, satellite and ground frames")
  System_Ext(tools, "Pattern sources", "TICRA GRASP, measurement ranges (GRD, CSV files)")
  System_Ext(github, "GitHub", "Hosts the sources; runs CI")
  Rel(developer, antenna, "Calls the API of", "Rust")
  Rel(antenna, rfr, "Uses dB and aperture formulas from", "Rust crate")
  Rel(antenna, geo, "Uses frames and vectors from", "Rust crate")
  Rel(developer, tools, "Exports pattern files from")
  Rel(developer, github, "Fetches the crates from", "Git over HTTPS")
```

Elements:
- **Application developer**: anyone building antenna or RF tooling in Rust.
- **antenna_r**: this crate.
- **rf-r**, **geometry_r**: sibling crates, git dependencies.
- **Pattern sources**: produce the GRD and CSV files the developer passes to the parsers.
- **GitHub**: source hosting and CI, build time only.
