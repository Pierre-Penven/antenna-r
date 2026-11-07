//! Simple tabular pattern input: a θ/φ gain table and a rotationally
//! symmetric 1-D cut. Both produce a [`GridPattern`] on a θ/φ grid.
//!
//! Format: one record per line, values separated by commas, semicolons or
//! whitespace; `#` starts a comment line; the first non-comment line may be
//! a text header. Angles are in degrees, gains in dBi.

use crate::error::{AntennaError, Result};
use crate::grd::{GridAxis, GridPattern, GridType};

/// Relative tolerance when checking that samples lie on a regular grid.
const REGULAR_TOLERANCE: f64 = 1e-6;

/// Parses `θ_deg, φ_deg, gain_dBi` records covering a full regular θ/φ
/// grid (any order).
///
/// ```
/// use antenna_r::{Direction, GainPattern, grd::table::parse_theta_phi_table};
/// let csv = "theta,phi,gain\n0,0,10\n0,90,10\n10,0,4\n10,90,6\n";
/// let p = parse_theta_phi_table(csv).unwrap();
/// let g = p.gain_dbi(Direction::from_theta_phi(5f64.to_radians(), 45f64.to_radians()));
/// assert!((g - 7.5).abs() < 1e-9);
/// ```
///
/// # Errors
/// [`AntennaError::Parse`] for malformed records; [`AntennaError::InvalidParameter`]
/// when the samples do not form a complete regular grid.
pub fn parse_theta_phi_table(text: &str) -> Result<GridPattern> {
    let records = records::<3>(text)?;
    let thetas = regular_axis("theta", records.iter().map(|r| r[0]))?;
    let phis = regular_axis("phi", records.iter().map(|r| r[1]))?;
    let mut gains = vec![f64::NAN; thetas.count * phis.count];
    for r in &records {
        let i = index_on(&phis, r[1]);
        let j = index_on(&thetas, r[0]);
        gains[j * phis.count + i] = r[2];
    }
    if gains.iter().any(|g| g.is_nan()) {
        return Err(crate::error::invalid(
            "table",
            format!(
                "{} records do not cover the {} × {} grid",
                records.len(),
                thetas.count,
                phis.count
            ),
        ));
    }
    GridPattern::new(GridType::ThetaPhi, phis, thetas, gains)
}

/// Parses a rotationally symmetric cut of `θ_deg, gain_dBi` records on a
/// regular θ grid (any order).
///
/// # Errors
/// As [`parse_theta_phi_table`].
pub fn parse_symmetric_cut(text: &str) -> Result<GridPattern> {
    let records = records::<2>(text)?;
    let thetas = regular_axis("theta", records.iter().map(|r| r[0]))?;
    let mut gains = vec![f64::NAN; thetas.count];
    for r in &records {
        gains[index_on(&thetas, r[0])] = r[1];
    }
    let phi = GridAxis::new(0.0, 0.0, 1)?;
    GridPattern::new(GridType::ThetaPhi, phi, thetas, gains)
}

fn records<const N: usize>(text: &str) -> Result<Vec<[f64; N]>> {
    let mut out = Vec::new();
    let mut first = true;
    for (n, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parsed = parse_record::<N>(trimmed);
        match parsed {
            Some(r) => out.push(r),
            None if first => {} // text header
            None => {
                return Err(AntennaError::Parse {
                    line: n + 1,
                    reason: format!("expected {N} finite numbers"),
                });
            }
        }
        first = false;
    }
    if out.is_empty() {
        return Err(AntennaError::Parse {
            line: text.lines().count() + 1,
            reason: "no records".into(),
        });
    }
    Ok(out)
}

fn parse_record<const N: usize>(line: &str) -> Option<[f64; N]> {
    let mut out = [0.0; N];
    let mut fields = line
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|f| !f.is_empty());
    for o in &mut out {
        *o = fields.next()?.parse().ok().filter(|v: &f64| v.is_finite())?;
    }
    fields.next().is_none().then_some(out)
}

/// Builds the axis of the distinct values, which must be regularly spaced.
fn regular_axis(name: &'static str, values: impl Iterator<Item = f64>) -> Result<GridAxis> {
    let mut v: Vec<f64> = values.collect();
    v.sort_by(f64::total_cmp);
    v.dedup();
    let axis = GridAxis::new(v[0], v[v.len() - 1], v.len())?;
    let step = axis.step();
    let tol = REGULAR_TOLERANCE * step.max(1.0);
    if v.iter().enumerate().any(|(i, x)| (x - axis.value(i)).abs() > tol) {
        return Err(crate::error::invalid(name, "values are not regularly spaced"));
    }
    Ok(axis)
}

fn index_on(axis: &GridAxis, value: f64) -> usize {
    if axis.count == 1 {
        0
    } else {
        ((value - axis.start) / axis.step()).round() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::direction::Direction;
    use crate::math::assert_close;
    use crate::pattern::GainPattern;

    #[test]
    fn table_accepts_any_order_separators_and_comments() {
        let text = "# measured\n10;90;6\n0 0 10\n\n0,90,10\n10,0,4\n";
        let p = parse_theta_phi_table(text).unwrap();
        assert_eq!(p.peak_gain_dbi(), 10.0);
        let d = Direction::from_theta_phi(10f64.to_radians(), 0.0);
        assert_close(p.gain_dbi(d), 4.0, 1e-9);
    }

    #[test]
    fn table_rejects_incomplete_or_irregular_grids() {
        assert!(matches!(
            parse_theta_phi_table("0,0,1\n0,90,1\n10,0,1\n"),
            Err(AntennaError::InvalidParameter { name: "table", .. })
        ));
        assert!(matches!(
            parse_theta_phi_table("0,0,1\n1,0,1\n3,0,1\n"),
            Err(AntennaError::InvalidParameter { name: "theta", .. })
        ));
    }

    #[test]
    fn table_reports_bad_lines() {
        assert!(matches!(
            parse_theta_phi_table("0,0,1\n0,x,1\n"),
            Err(AntennaError::Parse { line: 2, .. })
        ));
        assert!(matches!(
            parse_theta_phi_table("0,0,1,5\n"),
            Err(AntennaError::Parse { .. })
        ));
        assert!(matches!(
            parse_theta_phi_table("header only\n# c\n"),
            Err(AntennaError::Parse { line: 3, .. })
        ));
        assert!(matches!(
            parse_theta_phi_table(""),
            Err(AntennaError::Parse { line: 1, .. })
        ));
        assert!(matches!(
            parse_symmetric_cut("0,inf\n"),
            Err(AntennaError::Parse { .. })
        ));
    }

    #[test]
    fn symmetric_cut_ignores_azimuth() {
        let p = parse_symmetric_cut("theta gain\n0 20\n2 14\n4 0\n").unwrap();
        for phi in [0.0, 2.0, 5.0] {
            let d = Direction::from_theta_phi(1f64.to_radians(), phi);
            assert_close(p.gain_dbi(d), 17.0, 1e-9);
        }
        // Single-sample cut: a constant pattern at θ = 0 only.
        let single = parse_symmetric_cut("0 5\n").unwrap();
        assert_eq!(single.gain_dbi(Direction::BORESIGHT), 5.0);
    }
}
