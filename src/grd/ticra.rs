//! Reader and writer for TICRA GRASP far-field grid files (`.grd`, ASCII).
//!
//! Record layout (TICRA, "Field data in rectangular grid"):
//!
//! ```text
//! TEXT …                          header lines, until a line starting "++++"
//! KTYPE                           must be 1
//! NSET ICOMP NCOMP IGRID
//! IX IY                           NSET lines: beam centre offsets (cells)
//! per set:
//!   XS YS XE YE                   grid limits
//!   NX NY KLIMIT
//!   per row J = 1..NY:
//!     [IS IN]                     only when KLIMIT = 1
//!     F1re F1im F2re F2im [F3re F3im]   one line per point I = IS..IS+IN−1
//! ```
//!
//! Grid points: `X = XCEN + XS + DX (I − 1)`, `DX = (XE − XS)/(NX − 1)`,
//! `XCEN = DX · IX` (same for Y). Angles are in degrees.
//!
//! Field normalisation follows GRASP: `|F1|² + |F2|²` is the linear power
//! gain when the source power is normalised to 4π W (GRASP's usual
//! convention, assumed here). Use [`GridPattern::with_offset_db`] or
//! [`GridPattern::normalised_to_directivity`] when a file is normalised
//! differently.
//!
//! The input is treated as untrusted: every count is checked, the total
//! number of points is capped by [`MAX_GRID_POINTS`], and malformed input
//! returns [`AntennaError::Parse`] with the line number, never a panic.

use std::fmt::Write as _;

use crate::error::{AntennaError, Result};
use crate::grd::{GridAxis, GridPattern, GridType};
use crate::math::Complex;
use crate::pattern::GainPattern;

/// Largest total number of grid points accepted in one file (all sets).
pub const MAX_GRID_POINTS: usize = 1 << 24;
/// Largest number of field sets accepted in one file.
pub const MAX_SETS: usize = 4096;

/// Field components stored in a file (`ICOMP`). Only the codes whose total
/// power can be computed are supported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldComponents {
    /// `ICOMP = 1`: linear `Eθ`, `Eφ`.
    LinearThetaPhi,
    /// `ICOMP = 2`: right- and left-hand circular.
    Circular,
    /// `ICOMP = 3`: linear co and cross (Ludwig's third definition).
    Ludwig3,
    /// `ICOMP = 4`: major and minor axes of the polarisation ellipse
    /// (real parts, linear scale).
    PolarisationEllipse,
}

impl FieldComponents {
    /// From the TICRA `ICOMP` code.
    ///
    /// # Errors
    /// [`AntennaError::UnsupportedFormat`] for XPD ratios (5–8), total
    /// power with axial ratio (9) and unknown codes.
    pub fn from_icomp(code: i64) -> Result<Self> {
        match code {
            1 => Ok(Self::LinearThetaPhi),
            2 => Ok(Self::Circular),
            3 => Ok(Self::Ludwig3),
            4 => Ok(Self::PolarisationEllipse),
            _ => Err(AntennaError::UnsupportedFormat(format!("ICOMP = {code}"))),
        }
    }

    /// The TICRA `ICOMP` code.
    #[must_use]
    pub fn icomp(self) -> i64 {
        match self {
            Self::LinearThetaPhi => 1,
            Self::Circular => 2,
            Self::Ludwig3 => 3,
            Self::PolarisationEllipse => 4,
        }
    }

    /// Total linear power `|F1|² + |F2|²` of one point.
    #[must_use]
    pub fn total_power(self, f1: Complex, f2: Complex) -> f64 {
        match self {
            Self::PolarisationEllipse => f1.re.mul_add(f1.re, f2.re * f2.re),
            _ => f1.norm_sqr() + f2.norm_sqr(),
        }
    }
}

/// One field set (beam) of a grid file.
#[derive(Debug, Clone, PartialEq)]
pub struct GrdSet {
    /// X axis, beam centre offset included.
    pub x: GridAxis,
    /// Y axis, beam centre offset included.
    pub y: GridAxis,
    /// First field component per point, X fastest; `NaN` where absent.
    pub f1: Vec<Complex>,
    /// Second field component per point, X fastest; `NaN` where absent.
    pub f2: Vec<Complex>,
}

/// A parsed TICRA grid file.
#[derive(Debug, Clone, PartialEq)]
pub struct GrdFile {
    /// Identification text lines before `++++`.
    pub header: Vec<String>,
    /// Meaning of the two field components.
    pub components: FieldComponents,
    /// Grid coordinates.
    pub grid: GridType,
    /// Field sets (beams).
    pub sets: Vec<GrdSet>,
}

const MISSING: Complex = Complex::new(f64::NAN, f64::NAN);

impl GrdFile {
    /// Parses an ASCII grid file.
    ///
    /// # Errors
    /// [`AntennaError::Parse`] for malformed or truncated input or counts
    /// beyond the limits; [`AntennaError::UnsupportedFormat`] for `KTYPE ≠ 1`
    /// or unsupported `ICOMP` / `IGRID` codes.
    pub fn parse(text: &str) -> Result<Self> {
        let mut lines = Lines::new(text);
        let header = lines.header()?;
        let mut ktype = [0];
        let line = lines.integers(&mut ktype)?;
        if ktype[0] != 1 {
            return Err(AntennaError::UnsupportedFormat(format!(
                "KTYPE = {} at line {line}",
                ktype[0]
            )));
        }
        let mut control = [0; 4];
        let line = lines.integers(&mut control)?;
        let [nset, icomp, ncomp, igrid] = control;
        let nset = count_in(line, "NSET", nset, 1, MAX_SETS)?;
        let components = FieldComponents::from_icomp(icomp)?;
        let ncomp = count_in(line, "NCOMP", ncomp, 2, 3)?;
        let grid = GridType::from_igrid(igrid)?;
        let mut centres = Vec::with_capacity(nset);
        for _ in 0..nset {
            let mut centre = [0; 2];
            lines.integers(&mut centre)?;
            centres.push(centre);
        }
        let mut sets = Vec::with_capacity(nset);
        let mut total_points = 0usize;
        for centre in centres {
            let set = lines.set(centre, ncomp, &mut total_points)?;
            sets.push(set);
        }
        Ok(Self {
            header,
            components,
            grid,
            sets,
        })
    }

    /// Writes the file in TICRA ASCII format (`KLIMIT = 0`, `NCOMP = 2`,
    /// beam centres 0). Missing points are written as zero fields.
    #[must_use]
    pub fn to_grd_string(&self) -> String {
        let mut s = String::new();
        for h in &self.header {
            let _ = writeln!(s, "{h}");
        }
        let _ = writeln!(s, "++++");
        let _ = writeln!(s, "1");
        let _ = writeln!(
            s,
            "{} {} 2 {}",
            self.sets.len(),
            self.components.icomp(),
            self.grid.igrid()
        );
        for _ in &self.sets {
            let _ = writeln!(s, "0 0");
        }
        for set in &self.sets {
            let _ = writeln!(s, "{:E} {:E} {:E} {:E}", set.x.start, set.y.start, set.x.end, set.y.end);
            let _ = writeln!(s, "{} {} 0", set.x.count, set.y.count);
            for (a, b) in set.f1.iter().zip(&set.f2) {
                let z = |c: &Complex| if c.re.is_nan() { Complex::default() } else { *c };
                let (a, b) = (z(a), z(b));
                let _ = writeln!(s, "{:E} {:E} {:E} {:E}", a.re, a.im, b.re, b.im);
            }
        }
        s
    }

    /// Samples any pattern on a grid and stores it as a single-set file with
    /// `F1 = √G` (co-polar, Ludwig-3) and `F2 = 0`. Points without a
    /// direction (outside the visible u/v disc) are zero fields.
    #[must_use]
    pub fn from_pattern<P: GainPattern + ?Sized>(pattern: &P, grid: GridType, x: GridAxis, y: GridAxis) -> Self {
        let mut f1 = Vec::with_capacity(x.count * y.count);
        for j in 0..y.count {
            for i in 0..x.count {
                let amp = grid
                    .direction(x.value(i), y.value(j))
                    .map_or(0.0, |d| rf_r::db::amplitude_from_db(pattern.gain_dbi(d)));
                f1.push(Complex::new(amp, 0.0));
            }
        }
        let f2 = vec![Complex::default(); f1.len()];
        Self {
            header: vec!["antenna_r grid export".to_string()],
            components: FieldComponents::Ludwig3,
            grid,
            sets: vec![GrdSet { x, y, f1, f2 }],
        }
    }

    /// Gain pattern of set `index` (0-based), gain = `10 log10(|F1|² + |F2|²)`.
    ///
    /// # Errors
    /// [`AntennaError::OutOfRange`] for a missing set;
    /// [`AntennaError::InvalidParameter`] when the set has no point.
    pub fn pattern(&self, index: usize) -> Result<GridPattern> {
        let set = self
            .sets
            .get(index)
            .ok_or_else(|| AntennaError::OutOfRange(format!("set {index} of {}", self.sets.len())))?;
        let gains = set
            .f1
            .iter()
            .zip(&set.f2)
            .map(|(a, b)| rf_r::db::to_db(self.components.total_power(*a, *b)))
            .collect();
        GridPattern::new(self.grid, set.x, set.y, gains)
    }
}

fn too_few(line: usize, expected: usize) -> AntennaError {
    AntennaError::Parse {
        line,
        reason: format!("expected {expected} values, found fewer"),
    }
}

/// Validates an integer count within `[lo, hi]`.
fn count_in(line: usize, name: &str, value: i64, lo: usize, hi: usize) -> Result<usize> {
    usize::try_from(value)
        .ok()
        .filter(|v| (lo..=hi).contains(v))
        .ok_or_else(|| AntennaError::Parse {
            line,
            reason: format!("{name} = {value} outside [{lo}, {hi}]"),
        })
}

/// Parses a Fortran real: accepts `D`/`d` exponents and a missing exponent
/// letter (`1.0-05`).
fn parse_real(token: &str) -> Option<f64> {
    let t = token.replace(['D', 'd'], "E");
    if let Ok(v) = t.parse::<f64>() {
        return Some(v);
    }
    let split = t
        .char_indices()
        .skip(1)
        .find(|&(i, c)| (c == '+' || c == '-') && !matches!(t.as_bytes()[i - 1], b'E' | b'e'))?;
    format!("{}E{}", &t[..split.0], &t[split.0..]).parse().ok()
}

/// Line cursor over the non-blank lines of the input, with 1-based numbers.
struct Lines<'a> {
    inner: std::iter::Enumerate<std::str::Lines<'a>>,
    last_line: usize,
}

impl<'a> Lines<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            inner: text.lines().enumerate(),
            last_line: 0,
        }
    }

    fn next_line(&mut self) -> Result<(usize, &'a str)> {
        for (i, l) in self.inner.by_ref() {
            self.last_line = i + 1;
            if !l.trim().is_empty() {
                return Ok((i + 1, l));
            }
        }
        Err(AntennaError::Parse {
            line: self.last_line + 1,
            reason: "unexpected end of file".into(),
        })
    }

    fn header(&mut self) -> Result<Vec<String>> {
        let mut header = Vec::new();
        loop {
            let (_, l) = self.next_line()?;
            if l.starts_with("++++") {
                return Ok(header);
            }
            header.push(l.to_string());
        }
    }

    /// Parses the first `out.len()` fields of the next line as integers.
    fn integers(&mut self, out: &mut [i64]) -> Result<usize> {
        let (line, l) = self.next_line()?;
        let mut fields = l.split_whitespace();
        let expected = out.len();
        for o in out.iter_mut() {
            let t = fields.next().ok_or_else(|| too_few(line, expected))?;
            *o = t.parse().map_err(|_| AntennaError::Parse {
                line,
                reason: format!("expected an integer, found `{t}`"),
            })?;
        }
        Ok(line)
    }

    /// Parses the first `out.len()` fields of the next line as finite reals.
    fn reals(&mut self, out: &mut [f64]) -> Result<usize> {
        let (line, l) = self.next_line()?;
        let mut fields = l.split_whitespace();
        let expected = out.len();
        for o in out.iter_mut() {
            let t = fields.next().ok_or_else(|| too_few(line, expected))?;
            *o = parse_real(t)
                .filter(|v| v.is_finite())
                .ok_or_else(|| AntennaError::Parse {
                    line,
                    reason: format!("expected a finite real, found `{t}`"),
                })?;
        }
        Ok(line)
    }

    fn set(&mut self, centre: [i64; 2], ncomp: usize, total_points: &mut usize) -> Result<GrdSet> {
        let mut limits = [0.0; 4];
        let limits_line = self.reals(&mut limits)?;
        let [xs, ys, xe, ye] = limits;
        let mut counts = [0; 3];
        let line = self.integers(&mut counts)?;
        let [nx, ny, klimit] = counts;
        let nx = count_in(line, "NX", nx, 1, MAX_GRID_POINTS)?;
        let ny = count_in(line, "NY", ny, 1, MAX_GRID_POINTS)?;
        let points = nx
            .checked_mul(ny)
            .and_then(|p| p.checked_add(*total_points))
            .filter(|p| *p <= MAX_GRID_POINTS)
            .ok_or_else(|| AntennaError::Parse {
                line,
                reason: format!("grid too large (more than {MAX_GRID_POINTS} points)"),
            })?;
        *total_points = points;
        let limited = match klimit {
            0 => false,
            1 => true,
            _ => {
                return Err(AntennaError::Parse {
                    line,
                    reason: format!("KLIMIT = {klimit}, expected 0 or 1"),
                });
            }
        };
        let x = axis_with_centre(limits_line, xs, xe, nx, centre[0])?;
        let y = axis_with_centre(limits_line, ys, ye, ny, centre[1])?;
        let mut f1 = vec![MISSING; nx * ny];
        let mut f2 = vec![MISSING; nx * ny];
        for j in 0..ny {
            let (first, count) = if limited { self.row_limits(nx)? } else { (0, nx) };
            for i in first..first + count {
                let (f, g) = self.point(ncomp)?;
                f1[j * nx + i] = f;
                f2[j * nx + i] = g;
            }
        }
        Ok(GrdSet { x, y, f1, f2 })
    }

    /// Reads `IS IN` and returns the 0-based first column and the count.
    fn row_limits(&mut self, nx: usize) -> Result<(usize, usize)> {
        let mut limits = [0; 2];
        let line = self.integers(&mut limits)?;
        let [is, n] = limits;
        let bad = || AntennaError::Parse {
            line,
            reason: format!("row limits IS = {is}, IN = {n} outside 1..={nx}"),
        };
        let first = usize::try_from(is).ok().filter(|v| *v >= 1).ok_or_else(bad)? - 1;
        let count = usize::try_from(n).map_err(|_| bad())?;
        if first + count > nx {
            return Err(bad());
        }
        Ok((first, count))
    }

    fn point(&mut self, ncomp: usize) -> Result<(Complex, Complex)> {
        // The optional third (radial) component of near fields is ignored.
        let mut v = [0.0; 6];
        self.reals(&mut v[..2 * ncomp])?;
        Ok((Complex::new(v[0], v[1]), Complex::new(v[2], v[3])))
    }
}

/// Axis `XCEN + XS … XCEN + XE` with `XCEN = DX · IX`.
fn axis_with_centre(line: usize, start: f64, end: f64, count: usize, centre: i64) -> Result<GridAxis> {
    let step = if count > 1 {
        (end - start) / (count - 1) as f64
    } else {
        0.0
    };
    let offset = step * centre as f64;
    GridAxis::new(offset + start, offset + end, count).map_err(|e| AntennaError::Parse {
        line,
        reason: e.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::direction::Direction;
    use crate::math::assert_close;
    use crate::pattern::CosineQ;

    const SIMPLE: &str = "GRASP test file\nsecond header line\n++++\n1\n\n1 1 2 7\n0 0\n\
0.0 0.0 90.0 10.0\n2 2 0\n1.0 0.0 0.0 0.0\n0.0 1.0 0.0 0.0\n2.0 0.0 0.0 0.0\n1.0D+00 0.0 1.0 0.0\n";

    #[test]
    fn parses_header_counts_and_values() {
        let f = GrdFile::parse(SIMPLE).unwrap();
        assert_eq!(f.header, vec!["GRASP test file", "second header line"]);
        assert_eq!(f.components, FieldComponents::LinearThetaPhi);
        assert_eq!(f.grid, GridType::ThetaPhi);
        let s = &f.sets[0];
        assert_eq!(s.x, GridAxis::new(0.0, 90.0, 2).unwrap());
        assert_eq!(s.y, GridAxis::new(0.0, 10.0, 2).unwrap());
        assert_eq!(s.f1[1], Complex::new(0.0, 1.0));
        assert_eq!(s.f2[3], Complex::new(1.0, 0.0));
        let p = f.pattern(0).unwrap();
        assert_close(p.values_dbi()[0], 0.0, 1e-12);
        assert_close(p.values_dbi()[2], 10.0 * 4f64.log10(), 1e-12);
        assert_close(p.values_dbi()[3], 10.0 * 2f64.log10(), 1e-12);
    }

    #[test]
    fn beam_centre_offsets_shift_the_axes() {
        let text = "h\n++++\n1\n1 3 2 1\n2 -1\n-0.1 -0.1 0.1 0.1\n3 3 0\n".to_string() + &"1 0 0 0\n".repeat(9);
        let f = GrdFile::parse(&text).unwrap();
        assert_eq!(f.components, FieldComponents::Ludwig3);
        let s = &f.sets[0];
        assert_close(s.x.start, 0.1, 1e-15);
        assert_close(s.x.end, 0.3, 1e-15);
        assert_close(s.y.start, -0.2, 1e-15);
    }

    #[test]
    fn single_column_file_is_a_symmetric_cut() {
        let text = "h\n++++\n1\n1 1 2 7\n5 0\n0 0 0 10\n1 2 0\n2 0 0 0\n1 0 0 0\n";
        let f = GrdFile::parse(text).unwrap();
        // A single column has no spacing, so the centre offset is ignored.
        assert_eq!(f.sets[0].x, GridAxis::new(0.0, 0.0, 1).unwrap());
        let p = f.pattern(0).unwrap();
        let d = Direction::from_theta_phi(5f64.to_radians(), 2.0);
        assert_close(p.gain_dbi(d), 0.5 * (10.0 * 4f64.log10()), 1e-9);
    }

    #[test]
    fn klimit_rows_leave_missing_points() {
        let text = "h\n++++\n1\n1 2 3 1\n0 0\n-0.1 -0.1 0.1 0.1\n3 2 1\n2 1\n1 0 0 0 9 9\n1 3\n1 0 0 0 0 0\n2 0 0 0 0 0\n3 0 0 0 0 0\n";
        let f = GrdFile::parse(text).unwrap();
        assert_eq!(f.components, FieldComponents::Circular);
        let s = &f.sets[0];
        assert!(s.f1[0].re.is_nan());
        assert_eq!(s.f1[1], Complex::new(1.0, 0.0));
        assert!(s.f1[2].re.is_nan());
        assert_eq!(s.f1[5], Complex::new(3.0, 0.0));
        let p = f.pattern(0).unwrap();
        assert!(p.values_dbi()[0].is_nan());
    }

    #[test]
    fn klimit_row_with_no_points() {
        let text = "h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 1\n1 0\n1 2\n1 0 0 0\n1 0 0 0\n";
        let f = GrdFile::parse(text).unwrap();
        assert!(f.sets[0].f1[0].re.is_nan());
        assert_eq!(f.sets[0].f1[3], Complex::new(1.0, 0.0));
    }

    #[test]
    fn polarisation_ellipse_power_uses_real_parts() {
        let c = FieldComponents::PolarisationEllipse;
        assert_close(
            c.total_power(Complex::new(3.0, 9.0), Complex::new(4.0, 9.0)),
            25.0,
            1e-12,
        );
        assert_close(
            FieldComponents::Circular.total_power(Complex::new(3.0, 4.0), Complex::new(0.0, 1.0)),
            26.0,
            1e-12,
        );
        for code in 1..=4 {
            assert_eq!(FieldComponents::from_icomp(code).unwrap().icomp(), code);
        }
    }

    #[test]
    fn unsupported_codes_are_reported() {
        let bad_icomp = SIMPLE.replace("1 1 2 7", "1 5 2 7");
        assert!(matches!(
            GrdFile::parse(&bad_icomp),
            Err(AntennaError::UnsupportedFormat(_))
        ));
        let bad_igrid = SIMPLE.replace("1 1 2 7", "1 1 2 3");
        assert!(matches!(
            GrdFile::parse(&bad_igrid),
            Err(AntennaError::UnsupportedFormat(_))
        ));
        let bad_ktype = SIMPLE.replace("++++\n1\n", "++++\n2\n");
        assert!(matches!(
            GrdFile::parse(&bad_ktype),
            Err(AntennaError::UnsupportedFormat(_))
        ));
    }

    #[test]
    fn malformed_input_reports_the_line() {
        let cases = [
            ("no separator\n", 2),
            ("h\n++++\nx\n", 3),
            ("h\n++++\n1\n1 1 2\n", 4),
            ("h\n++++\n1\n0 1 2 1\n", 4),
            ("h\n++++\n1\n1 1 4 1\n", 4),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1\n", 6),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 inf\n", 6),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 2\n", 7),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n0 2 0\n", 7),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 0 1\n2 2 0\n", 6),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 1\n0 2\n", 8),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 1\n2 2\n", 8),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 1\n1 -1\n", 8),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 0\n1 0 0\n", 8),
            ("h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n2 2 0\n1 0 0 0\n", 9),
        ];
        for (text, line) in cases {
            let err = GrdFile::parse(text).unwrap_err();
            assert!(
                matches!(err, AntennaError::Parse { line: l, .. } if l == line),
                "{text:?}: {err:?}"
            );
        }
    }

    #[test]
    fn oversized_grids_are_rejected_before_allocation() {
        let text = "h\n++++\n1\n1 1 2 1\n0 0\n0 0 1 1\n100000 100000 0\n";
        let err = GrdFile::parse(text).unwrap_err();
        assert!(matches!(&err, AntennaError::Parse { line: 7, reason } if reason.contains("too large")));
        let too_many_sets = "h\n++++\n1\n5000 1 2 1\n";
        assert!(matches!(
            GrdFile::parse(too_many_sets),
            Err(AntennaError::Parse { line: 4, .. })
        ));
    }

    #[test]
    fn fortran_reals_are_accepted() {
        assert_eq!(parse_real("1.5D+01"), Some(15.0));
        assert_eq!(parse_real("2.5d-1"), Some(0.25));
        assert_eq!(parse_real("1.0-05"), Some(1e-5));
        assert_eq!(parse_real("-1.0+02"), Some(-100.0));
        assert_eq!(parse_real("1.0E-05"), Some(1e-5));
        assert_eq!(parse_real("abc"), None);
        assert_eq!(parse_real("1.0e+-5"), None);
    }

    #[test]
    fn missing_set_is_out_of_range() {
        let f = GrdFile::parse(SIMPLE).unwrap();
        assert!(matches!(f.pattern(1), Err(AntennaError::OutOfRange(_))));
    }

    #[test]
    fn writer_round_trips_through_the_parser() {
        let p = CosineQ::new(2.0).unwrap();
        let x = GridAxis::new(0.0, 350.0, 36).unwrap();
        let y = GridAxis::new(0.0, 180.0, 37).unwrap();
        let file = GrdFile::from_pattern(&p, GridType::ThetaPhi, x, y);
        let text = file.to_grd_string();
        let back = GrdFile::parse(&text).unwrap();
        assert_eq!(back, file);
        let grid = back.pattern(0).unwrap();
        for (theta, phi) in [(0.0, 0.0), (0.3, 1.0), (1.2, 4.0)] {
            let d = Direction::from_theta_phi(theta, phi);
            assert_close(grid.gain_dbi(Direction::from_theta_phi(theta, phi)), p.gain_dbi(d), 0.2);
        }
        // Exact at grid nodes.
        let node = Direction::from_theta_phi(30f64.to_radians(), 0.0);
        assert_close(grid.gain_dbi(node), p.gain_dbi(node), 1e-9);
    }

    #[test]
    fn writer_zeroes_missing_and_invisible_points() {
        let p = CosineQ::new(1.0).unwrap();
        let ax = GridAxis::new(-1.0, 1.0, 3).unwrap();
        let mut file = GrdFile::from_pattern(&p, GridType::Uv, ax, ax);
        assert_eq!(file.sets[0].f1[0], Complex::default());
        file.sets[0].f1[1] = MISSING;
        let back = GrdFile::parse(&file.to_grd_string()).unwrap();
        assert_eq!(back.sets[0].f1[1], Complex::default());
    }
}
