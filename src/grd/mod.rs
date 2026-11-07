//! Gain patterns sampled on a regular 2-D grid: TICRA GRASP `.grd` files
//! ([`ticra`]), CSV tables and 1-D cuts ([`table`]), and the
//! [`GridPattern`] that interpolates them.

pub mod table;
pub mod ticra;

use crate::analysis::mean_gain_db_piecewise;
use crate::direction::Direction;
use crate::error::{Result, invalid, require_finite};
use crate::pattern::GainPattern;

/// Lowest gain stored in a grid (dBi): zero fields become this value so
/// interpolation in dB stays finite.
pub const MIN_GRID_GAIN_DBI: f64 = -300.0;

/// Tolerance (in grid cells) for points that sit exactly on a grid edge.
const EDGE_EPS: f64 = 1e-9;

/// Full turn in degrees, the period of azimuth-like grid axes.
const FULL_TURN_DEG: f64 = 360.0;

/// Coordinates of a grid, with the TICRA GRASP `IGRID` definitions. Angles
/// are in degrees; `u`, `v` are direction cosines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridType {
    /// `IGRID = 1`: `(X, Y) = (u, v)`, front hemisphere.
    Uv,
    /// `IGRID = 4`: elevation over azimuth, `(X, Y) = (Az, El)`,
    /// `r̂ = (−sin Az cos El, sin El, cos Az cos El)`.
    ElOverAz,
    /// `IGRID = 5`: elevation and azimuth, `(X, Y) = (Az, El)`,
    /// `Az = −θ cos φ`, `El = θ sin φ`.
    ElAndAz,
    /// `IGRID = 6`: azimuth over elevation, `(X, Y) = (Az, El)`,
    /// `r̂ = (−sin Az, cos Az sin El, cos Az cos El)`.
    AzOverEl,
    /// `IGRID = 7`: `(X, Y) = (φ, θ)`.
    ThetaPhi,
}

impl GridType {
    /// From the TICRA `IGRID` code.
    ///
    /// # Errors
    /// [`AntennaError::UnsupportedFormat`](crate::AntennaError::UnsupportedFormat)
    /// for any other code.
    pub fn from_igrid(code: i64) -> Result<Self> {
        match code {
            1 => Ok(Self::Uv),
            4 => Ok(Self::ElOverAz),
            5 => Ok(Self::ElAndAz),
            6 => Ok(Self::AzOverEl),
            7 => Ok(Self::ThetaPhi),
            _ => Err(crate::AntennaError::UnsupportedFormat(format!("IGRID = {code}"))),
        }
    }

    /// The TICRA `IGRID` code.
    #[must_use]
    pub fn igrid(self) -> i64 {
        match self {
            Self::Uv => 1,
            Self::ElOverAz => 4,
            Self::ElAndAz => 5,
            Self::AzOverEl => 6,
            Self::ThetaPhi => 7,
        }
    }

    /// Grid coordinates `(X, Y)` of a direction. For [`GridType::Uv`] the
    /// back hemisphere has no coordinates (`None`).
    #[must_use]
    pub fn coordinates(self, dir: Direction) -> Option<(f64, f64)> {
        let deg = |(a, b): (f64, f64)| (a.to_degrees(), b.to_degrees());
        match self {
            Self::Uv => (dir.w() >= 0.0).then(|| (dir.u(), dir.v())),
            Self::ElOverAz => Some(deg(dir.to_el_over_az())),
            Self::ElAndAz => Some(deg(dir.to_el_and_az())),
            Self::AzOverEl => Some(deg(dir.to_az_over_el())),
            Self::ThetaPhi => Some((dir.phi_rad().to_degrees(), dir.theta_rad().to_degrees())),
        }
    }

    /// Direction of grid point `(x, y)`; `None` outside the visible region
    /// of a u/v grid.
    #[must_use]
    pub fn direction(self, x: f64, y: f64) -> Option<Direction> {
        match self {
            Self::Uv => Direction::from_uv(x, y).ok(),
            Self::ElOverAz => Some(Direction::from_el_over_az(x.to_radians(), y.to_radians())),
            Self::ElAndAz => Some(Direction::from_el_and_az(x.to_radians(), y.to_radians())),
            Self::AzOverEl => Some(Direction::from_az_over_el(x.to_radians(), y.to_radians())),
            Self::ThetaPhi => Some(Direction::from_theta_phi(y.to_radians(), x.to_radians())),
        }
    }

    /// Whether the X / Y coordinate is an angle that wraps every 360°.
    fn periodic_axes(self) -> (bool, bool) {
        match self {
            Self::ElOverAz | Self::ThetaPhi => (true, false),
            Self::AzOverEl => (false, true),
            Self::Uv | Self::ElAndAz => (false, false),
        }
    }
}

/// One regular grid axis: `count` points from `start` to `end` inclusive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridAxis {
    /// First coordinate.
    pub start: f64,
    /// Last coordinate (equal to `start` when `count == 1`).
    pub end: f64,
    /// Number of points (≥ 1).
    pub count: usize,
}

impl GridAxis {
    /// Validated axis.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for non-finite limits, `count == 0`, or `end ≤ start` with `count > 1`.
    pub fn new(start: f64, end: f64, count: usize) -> Result<Self> {
        require_finite("start", start)?;
        require_finite("end", end)?;
        if count == 0 {
            return Err(invalid("count", "an axis needs at least one point"));
        }
        if count > 1 && end <= start {
            return Err(invalid("end", format!("must exceed start ({start}), got {end}")));
        }
        Ok(Self { start, end, count })
    }

    /// Spacing between points (0 for a single point).
    #[must_use]
    pub fn step(&self) -> f64 {
        if self.count > 1 {
            (self.end - self.start) / (self.count - 1) as f64
        } else {
            0.0
        }
    }

    /// Coordinate of point `i`.
    #[must_use]
    pub fn value(&self, i: usize) -> f64 {
        self.start + self.step() * i as f64
    }

    /// Whether the axis covers a full turn (the last cell wraps to the first).
    fn covers_full_turn(&self) -> bool {
        self.count > 1 && self.end - self.start + self.step() >= FULL_TURN_DEG - 1e-9
    }

    /// Cell index and fraction of coordinate `x`, or `None` outside.
    /// `wrap` closes the last cell onto the first for full-turn axes.
    fn locate(&self, x: f64, wrap: bool, clamp: bool) -> Option<(usize, usize, f64)> {
        if self.count == 1 {
            return ((x - self.start).abs() <= EDGE_EPS || clamp).then_some((0, 0, 0.0));
        }
        let last = self.count - 1;
        let mut f = (x - self.start) / self.step();
        if wrap && f > last as f64 && f < self.count as f64 {
            return Some((last, 0, f - last as f64));
        }
        if clamp {
            f = f.clamp(0.0, last as f64);
        }
        if !(f >= -EDGE_EPS && f <= last as f64 + EDGE_EPS) {
            return None;
        }
        let i = (f.floor().max(0.0) as usize).min(last - 1);
        Some((i, i + 1, (f - i as f64).clamp(0.0, 1.0)))
    }
}

/// How values between grid points are interpolated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Interpolation {
    /// Bilinear on gain in dB (smooth near deep nulls). Default.
    #[default]
    BilinearDb,
    /// Bilinear on linear power gain.
    BilinearLinear,
}

/// What a grid returns for directions it does not cover.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum OutOfGrid {
    /// The lowest gain of the grid. Default.
    #[default]
    GridMinimum,
    /// A fixed gain in dBi.
    Floor(f64),
    /// The value at the nearest grid edge.
    Clamp,
    /// `NaN`, so the caller can detect it.
    Nan,
}

/// Gain pattern interpolated on a regular grid.
#[derive(Debug, Clone, PartialEq)]
pub struct GridPattern {
    grid: GridType,
    x: GridAxis,
    y: GridAxis,
    gain_dbi: Vec<f64>,
    min_gain_dbi: f64,
    peak_gain_dbi: f64,
    interpolation: Interpolation,
    out_of_grid: OutOfGrid,
    wrap_x: bool,
    wrap_y: bool,
}

impl GridPattern {
    /// Grid pattern from gains in dBi, X varying fastest:
    /// `gain_dbi[j * x.count + i]`. `NaN` marks a missing point; gains
    /// below [`MIN_GRID_GAIN_DBI`] (including `-inf`) are raised to it.
    ///
    /// ```
    /// use antenna_r::{Direction, GainPattern};
    /// use antenna_r::grd::{GridAxis, GridPattern, GridType};
    /// let x = GridAxis::new(-0.1, 0.1, 3).unwrap();
    /// let y = GridAxis::new(-0.1, 0.1, 3).unwrap();
    /// let g = vec![0.0, 3.0, 0.0, 3.0, 6.0, 3.0, 0.0, 3.0, 0.0];
    /// let p = GridPattern::new(GridType::Uv, x, y, g).unwrap();
    /// assert_eq!(p.gain_dbi(Direction::BORESIGHT), 6.0);
    /// assert!((p.gain_dbi(Direction::from_uv(0.05, 0.0).unwrap()) - 4.5).abs() < 1e-12);
    /// ```
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// when the value count is not `x.count · y.count`, a value is `+inf`,
    /// or no value is present.
    pub fn new(grid: GridType, x: GridAxis, y: GridAxis, mut gain_dbi: Vec<f64>) -> Result<Self> {
        if gain_dbi.len() != x.count * y.count {
            return Err(invalid(
                "gain_dbi",
                format!("expected {} values, got {}", x.count * y.count, gain_dbi.len()),
            ));
        }
        if gain_dbi.contains(&f64::INFINITY) {
            return Err(invalid("gain_dbi", "contains +inf"));
        }
        for g in &mut gain_dbi {
            if *g < MIN_GRID_GAIN_DBI {
                *g = MIN_GRID_GAIN_DBI;
            }
        }
        let present = gain_dbi.iter().copied().filter(|g| !g.is_nan());
        let (min, max) = present.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), g| (lo.min(g), hi.max(g)));
        if min > max {
            return Err(invalid("gain_dbi", "no grid value is present"));
        }
        let (px, py) = grid.periodic_axes();
        Ok(Self {
            grid,
            x,
            y,
            gain_dbi,
            min_gain_dbi: min,
            peak_gain_dbi: max,
            interpolation: Interpolation::default(),
            out_of_grid: OutOfGrid::default(),
            wrap_x: px,
            wrap_y: py,
        })
    }

    /// Selects the interpolation.
    #[must_use]
    pub fn with_interpolation(mut self, interpolation: Interpolation) -> Self {
        self.interpolation = interpolation;
        self
    }

    /// Selects the out-of-grid policy.
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a `Floor` that is not finite.
    pub fn with_out_of_grid(mut self, policy: OutOfGrid) -> Result<Self> {
        if let OutOfGrid::Floor(f) = policy {
            require_finite("floor", f)?;
        }
        self.out_of_grid = policy;
        Ok(self)
    }

    /// Adds `offset_db` to every value (e.g. a loss or a normalisation).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// for a non-finite offset.
    pub fn with_offset_db(mut self, offset_db: f64) -> Result<Self> {
        require_finite("offset_db", offset_db)?;
        for g in &mut self.gain_dbi {
            *g += offset_db;
        }
        self.min_gain_dbi += offset_db;
        self.peak_gain_dbi += offset_db;
        Ok(self)
    }

    /// Rescales the grid so that its power integrates to 0 dB over the
    /// sphere, i.e. the values become directivity (see
    /// [`mean_gain_db_piecewise`] for the integration parameters).
    ///
    /// # Errors
    /// [`AntennaError::InvalidParameter`](crate::AntennaError::InvalidParameter)
    /// if the integral is not finite.
    pub fn normalised_to_directivity(
        self,
        theta_edges_rad: &[f64],
        nodes_per_interval: usize,
        n_phi: usize,
    ) -> Result<Self> {
        let mean = mean_gain_db_piecewise(&self, theta_edges_rad, nodes_per_interval, n_phi);
        if !mean.is_finite() {
            return Err(invalid("gain_dbi", format!("power integral is not finite ({mean} dB)")));
        }
        self.with_offset_db(-mean)
    }

    /// Grid coordinates.
    #[must_use]
    pub fn grid_type(&self) -> GridType {
        self.grid
    }

    /// X and Y axes.
    #[must_use]
    pub fn axes(&self) -> (GridAxis, GridAxis) {
        (self.x, self.y)
    }

    /// Stored gains (dBi), X fastest.
    #[must_use]
    pub fn values_dbi(&self) -> &[f64] {
        &self.gain_dbi
    }

    /// Lowest stored gain (dBi).
    #[must_use]
    pub fn min_gain_dbi(&self) -> f64 {
        self.min_gain_dbi
    }

    fn outside(&self) -> f64 {
        match self.out_of_grid {
            OutOfGrid::GridMinimum | OutOfGrid::Clamp => self.min_gain_dbi,
            OutOfGrid::Floor(f) => f,
            OutOfGrid::Nan => f64::NAN,
        }
    }

    /// Wraps a periodic coordinate into `[start, start + 360)`.
    fn wrap(value: f64, axis: &GridAxis, periodic: bool) -> f64 {
        if periodic {
            axis.start + (value - axis.start).rem_euclid(FULL_TURN_DEG)
        } else {
            value
        }
    }

    fn sample(&self, x: f64, y: f64, clamp: bool) -> Option<f64> {
        let x = Self::wrap(x, &self.x, self.wrap_x);
        let y = Self::wrap(y, &self.y, self.wrap_y);
        let (i0, i1, tx) = self.x.locate(x, self.wrap_x && self.x.covers_full_turn(), clamp)?;
        let (j0, j1, ty) = self.y.locate(y, self.wrap_y && self.y.covers_full_turn(), clamp)?;
        let at = |i: usize, j: usize| self.gain_dbi[j * self.x.count + i];
        let corners = [at(i0, j0), at(i1, j0), at(i0, j1), at(i1, j1)];
        let weights = [(1.0 - tx) * (1.0 - ty), tx * (1.0 - ty), (1.0 - tx) * ty, tx * ty];
        // Zero-weight corners are skipped so a missing neighbour does not
        // poison a value that lies exactly on a grid line.
        let used = corners.iter().zip(weights).filter(|(_, w)| *w > 0.0);
        let value = match self.interpolation {
            Interpolation::BilinearDb => used.map(|(g, w)| g * w).sum::<f64>(),
            Interpolation::BilinearLinear => rf_r::db::to_db(used.map(|(g, w)| rf_r::db::from_db(*g) * w).sum::<f64>()),
        };
        (!value.is_nan()).then_some(value)
    }

    fn lookup(&self, dir: Direction, clamp: bool) -> Option<f64> {
        let (x, y) = self.grid.coordinates(dir)?;
        if self.grid == GridType::ThetaPhi {
            // Rotationally symmetric cut, or θ ∈ [−180°, 180°] grids that
            // store (θ, φ) as (−θ, φ + 180°).
            let x = if self.x.count == 1 { self.x.start } else { x };
            return self
                .sample(x, y, false)
                .or_else(|| self.sample(x + 180.0, -y, false))
                .or_else(|| if clamp { self.sample(x, y, true) } else { None });
        }
        self.sample(x, y, clamp)
    }
}

impl GainPattern for GridPattern {
    fn gain_dbi(&self, dir: Direction) -> f64 {
        if dir.w().is_nan() {
            return f64::NAN;
        }
        let clamp = self.out_of_grid == OutOfGrid::Clamp;
        self.lookup(dir, clamp).unwrap_or_else(|| self.outside())
    }
    fn peak_gain_dbi(&self) -> f64 {
        self.peak_gain_dbi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AntennaError;
    use crate::math::assert_close;
    use std::f64::consts::PI;

    fn uv_grid() -> GridPattern {
        let x = GridAxis::new(-0.1, 0.1, 3).unwrap();
        let y = GridAxis::new(-0.1, 0.1, 3).unwrap();
        GridPattern::new(
            GridType::Uv,
            x,
            y,
            vec![0.0, 3.0, 0.0, 3.0, 6.0, 3.0, 0.0, 3.0, f64::NAN],
        )
        .unwrap()
    }

    #[test]
    fn grid_type_codes_round_trip() {
        for g in [
            GridType::Uv,
            GridType::ElOverAz,
            GridType::ElAndAz,
            GridType::AzOverEl,
            GridType::ThetaPhi,
        ] {
            assert_eq!(GridType::from_igrid(g.igrid()).unwrap(), g);
        }
        assert!(matches!(
            GridType::from_igrid(2),
            Err(AntennaError::UnsupportedFormat(_))
        ));
    }

    #[test]
    fn grid_type_coordinates_round_trip() {
        let d = Direction::from_theta_phi(0.3, 2.0);
        for g in [
            GridType::Uv,
            GridType::ElOverAz,
            GridType::ElAndAz,
            GridType::AzOverEl,
            GridType::ThetaPhi,
        ] {
            let (x, y) = g.coordinates(d).unwrap();
            assert!(g.direction(x, y).unwrap().angle_to(d) < 1e-12, "{g:?}");
        }
        assert_eq!(GridType::Uv.coordinates(Direction::from_theta_phi(2.0, 0.0)), None);
        assert_eq!(GridType::Uv.direction(1.0, 1.0), None);
    }

    #[test]
    fn axis_validation_and_values() {
        let a = GridAxis::new(0.0, 10.0, 11).unwrap();
        assert_eq!(a.step(), 1.0);
        assert_eq!(a.value(3), 3.0);
        assert_eq!(GridAxis::new(5.0, 5.0, 1).unwrap().step(), 0.0);
        assert!(matches!(
            GridAxis::new(0.0, 1.0, 0),
            Err(AntennaError::InvalidParameter { name: "count", .. })
        ));
        assert!(matches!(
            GridAxis::new(1.0, 1.0, 2),
            Err(AntennaError::InvalidParameter { name: "end", .. })
        ));
        assert!(GridAxis::new(f64::NAN, 1.0, 2).is_err());
        assert!(GridAxis::new(0.0, f64::INFINITY, 2).is_err());
    }

    #[test]
    fn bilinear_interpolation_in_db_is_exact_at_nodes() {
        let p = uv_grid();
        assert_eq!(p.gain_dbi(Direction::BORESIGHT), 6.0);
        assert_eq!(p.peak_gain_dbi(), 6.0);
        assert_eq!(p.min_gain_dbi(), 0.0);
        assert_close(p.gain_dbi(Direction::from_uv(0.05, -0.05).unwrap()), 3.0, 1e-12);
        assert_close(p.gain_dbi(Direction::from_uv(-0.1, -0.1).unwrap()), 0.0, 1e-12);
    }

    #[test]
    fn bilinear_interpolation_in_linear_power() {
        let p = uv_grid().with_interpolation(Interpolation::BilinearLinear);
        let expected = 10.0 * f64::midpoint(10f64.powf(0.3), 10f64.powf(0.6)).log10();
        assert_close(p.gain_dbi(Direction::from_uv(0.05, 0.0).unwrap()), expected, 1e-12);
    }

    #[test]
    fn missing_points_and_outside_use_the_policy() {
        let p = uv_grid();
        let missing = Direction::from_uv(0.09, 0.09).unwrap();
        let outside = Direction::from_uv(0.5, 0.0).unwrap();
        let back = Direction::from_theta_phi(2.5, 0.0);
        assert_eq!(p.gain_dbi(missing), 0.0);
        assert_eq!(p.gain_dbi(outside), 0.0);
        assert_eq!(p.gain_dbi(back), 0.0);
        let p = p.with_out_of_grid(OutOfGrid::Floor(-20.0)).unwrap();
        assert_eq!(p.gain_dbi(outside), -20.0);
        let p = p.with_out_of_grid(OutOfGrid::Nan).unwrap();
        assert!(p.gain_dbi(outside).is_nan());
        let p = p.with_out_of_grid(OutOfGrid::Clamp).unwrap();
        assert_eq!(p.gain_dbi(Direction::from_uv(0.5, 0.0).unwrap()), 3.0);
        assert_eq!(p.gain_dbi(back), 0.0);
        assert!(p.gain_dbi(Direction::from_theta_phi(f64::NAN, 0.0)).is_nan());
        assert!(matches!(
            uv_grid().with_out_of_grid(OutOfGrid::Floor(f64::NAN)),
            Err(AntennaError::InvalidParameter { name: "floor", .. })
        ));
    }

    #[test]
    fn construction_validates_values() {
        let x = GridAxis::new(0.0, 1.0, 2).unwrap();
        assert!(matches!(
            GridPattern::new(GridType::Uv, x, x, vec![0.0; 3]),
            Err(AntennaError::InvalidParameter { name: "gain_dbi", .. })
        ));
        assert!(GridPattern::new(GridType::Uv, x, x, vec![f64::INFINITY; 4]).is_err());
        assert!(GridPattern::new(GridType::Uv, x, x, vec![f64::NAN; 4]).is_err());
        let p = GridPattern::new(GridType::Uv, x, x, vec![f64::NEG_INFINITY, 0.0, 0.0, 0.0]).unwrap();
        assert_eq!(p.values_dbi()[0], MIN_GRID_GAIN_DBI);
        assert_eq!(p.grid_type(), GridType::Uv);
        assert_eq!(p.axes(), (x, x));
    }

    #[test]
    fn offset_shifts_all_values() {
        let p = uv_grid().with_offset_db(-1.5).unwrap();
        assert_eq!(p.peak_gain_dbi(), 4.5);
        assert_eq!(p.min_gain_dbi(), -1.5);
        assert_eq!(p.gain_dbi(Direction::BORESIGHT), 4.5);
        assert!(uv_grid().with_offset_db(f64::NAN).is_err());
    }

    fn theta_phi_isotropic(phi_end: f64, n_phi: usize, value: f64) -> GridPattern {
        let x = GridAxis::new(0.0, phi_end, n_phi).unwrap();
        let y = GridAxis::new(0.0, 180.0, 19).unwrap();
        GridPattern::new(GridType::ThetaPhi, x, y, vec![value; n_phi * 19]).unwrap()
    }

    #[test]
    fn directivity_normalisation_brings_the_mean_to_zero() {
        let p = theta_phi_isotropic(350.0, 36, 7.0)
            .normalised_to_directivity(&[0.0, PI], 16, 8)
            .unwrap();
        assert_close(p.peak_gain_dbi(), 0.0, 1e-9);
        let x = GridAxis::new(0.0, 1.0, 2).unwrap();
        let empty = GridPattern::new(GridType::Uv, x, x, vec![MIN_GRID_GAIN_DBI; 4])
            .unwrap()
            .with_out_of_grid(OutOfGrid::Nan)
            .unwrap();
        assert!(empty.normalised_to_directivity(&[0.0, PI], 4, 4).is_err());
    }

    #[test]
    fn theta_phi_grid_wraps_the_last_azimuth_cell() {
        // φ = 0..350 by 10: φ = 355 interpolates between 350 and 0 (360).
        let x = GridAxis::new(0.0, 350.0, 36).unwrap();
        let y = GridAxis::new(0.0, 10.0, 2).unwrap();
        let mut g = vec![0.0; 72];
        g[35] = 10.0;
        g[71] = 10.0;
        let p = GridPattern::new(GridType::ThetaPhi, x, y, g).unwrap();
        let d = Direction::from_theta_phi(5f64.to_radians(), 355f64.to_radians());
        assert_close(p.gain_dbi(d), 5.0, 1e-9);
        let d = Direction::from_theta_phi(5f64.to_radians(), (-10f64).to_radians());
        assert_close(p.gain_dbi(d), 10.0, 1e-9);
    }

    #[test]
    fn theta_phi_grid_with_signed_theta_uses_the_mirror_point() {
        // θ ∈ [−90, 90], φ ∈ [0, 90]: (θ = 30°, φ = 200°) is (−30°, 20°).
        let x = GridAxis::new(0.0, 90.0, 10).unwrap();
        let y = GridAxis::new(-90.0, 90.0, 19).unwrap();
        let g: Vec<f64> = (0..19)
            .flat_map(|j| (0..10).map(move |i| f64::from(j * 100 + i)))
            .collect();
        let p = GridPattern::new(GridType::ThetaPhi, x, y, g).unwrap();
        let d = Direction::from_theta_phi(30f64.to_radians(), 200f64.to_radians());
        // j for −30° is 6, i for 20° is 2.
        assert_close(p.gain_dbi(d), 602.0, 1e-9);
        let d = Direction::from_theta_phi(30f64.to_radians(), 20f64.to_radians());
        assert_close(p.gain_dbi(d), 1202.0, 1e-9);
        // Outside both representations (θ = 120°) → clamp to the nearest edge.
        let p = p.with_out_of_grid(OutOfGrid::Clamp).unwrap();
        let d = Direction::from_theta_phi(120f64.to_radians(), 20f64.to_radians());
        assert_close(p.gain_dbi(d), 1802.0, 1e-9);
    }

    #[test]
    fn single_column_theta_phi_grid_is_rotationally_symmetric() {
        let x = GridAxis::new(0.0, 0.0, 1).unwrap();
        let y = GridAxis::new(0.0, 90.0, 4).unwrap();
        let p = GridPattern::new(GridType::ThetaPhi, x, y, vec![10.0, 7.0, 1.0, -5.0]).unwrap();
        for phi in [0.0, 1.0, 4.0] {
            let d = Direction::from_theta_phi(15f64.to_radians(), phi);
            assert_close(p.gain_dbi(d), 8.5, 1e-9);
        }
    }

    #[test]
    fn az_over_el_grid_wraps_elevation() {
        let x = GridAxis::new(-10.0, 10.0, 3).unwrap();
        let y = GridAxis::new(-180.0, 170.0, 36).unwrap();
        let p = GridPattern::new(GridType::AzOverEl, x, y, vec![1.0; 108]).unwrap();
        let d = Direction::from_az_over_el(0.0, 175f64.to_radians());
        assert_close(p.gain_dbi(d), 1.0, 1e-12);
    }

    #[test]
    fn single_point_axis_matches_only_its_coordinate() {
        let x = GridAxis::new(0.0, 0.0, 1).unwrap();
        let y = GridAxis::new(-0.1, 0.1, 3).unwrap();
        let p = GridPattern::new(GridType::Uv, x, y, vec![1.0, 2.0, 3.0])
            .unwrap()
            .with_out_of_grid(OutOfGrid::Nan)
            .unwrap();
        assert_close(p.gain_dbi(Direction::from_uv(0.0, 0.05).unwrap()), 2.5, 1e-12);
        assert!(p.gain_dbi(Direction::from_uv(0.01, 0.0).unwrap()).is_nan());
        let p = p.with_out_of_grid(OutOfGrid::Clamp).unwrap();
        assert_close(p.gain_dbi(Direction::from_uv(0.01, 0.0).unwrap()), 2.0, 1e-12);
    }
}
