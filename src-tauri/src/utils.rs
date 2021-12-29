pub fn lerp(p1: &ndarray::Array1<f32>, p2: &ndarray::Array1<f32>, t: f32) -> ndarray::Array1<f32> {
  p1 + (p2 - p1) * t
}
/// A struct that stores x and y
#[derive(Debug, Clone, Copy)]
pub struct Point<T> {
  pub x: T,
  pub y: T,
}
/// linear interpolation between 4 points
/// point: the point to interpolate to
///
/// origin: the point to be interpolated from that is closest to the world orignin (0,0)
///
/// vals: the values of the 4 points. Order should be `[(origin.x, origin.y),(origin.x, origin.y + 1),
/// (origin.x + 1, origin.y), (origin.x + 1, origin.y + 1)]`
///
///(0,0)
///    -------
///   | 0 | 2 |
///   ---------
///   | 1 | 3 |
///    -------
///
/// returns: the interpolated value
/// Idea from: https://blog.demofox.org/2015/04/30/bilinear-filtering-bilinear-interpolation/
pub fn square_lerp(point: Point<f32>, origin: Point<f32>, vals: ndarray::Array1<&f32>) -> f32 {
  let x_top = vals[1] + (point.x - origin.x) * (vals[3] - vals[1]);
  let x_bottom = vals[0] + (point.x - origin.x) * (vals[2] - vals[0]);
  x_bottom + (x_top - x_bottom) * (point.y - origin.y)
}

/// Color struct that stores red, green, and blue values
#[derive(Debug, Clone, Copy)]
pub struct Color<T> {
  pub r: T,
  pub g: T,
  pub b: T,
}
pub fn color_to_ndarray(color: Color<f32>) -> ndarray::Array1<f32> {
  ndarray::Array1::from(vec![color.r, color.g, color.b])
}
/// Linearly interpolate a color from a float between -1 and 1.
/// The if the value is outside of this range, it will be clamped to the closest value.
///
/// if val is -1, the color will be fully c1
/// .If the val is 1, the color will be fully c2
/// . If the val is 0, the color will be fully black
pub fn lerp_color(val: f32, c1: Color<f32>, c2: Color<f32>) -> Color<f32> {
  let val = val.clamp(-1.0, 1.0);
  let black = Color {
    r: 0.,
    g: 0.,
    b: 0.,
  };
  if val < 0. {
    let interpolated_color = lerp(&color_to_ndarray(black), &color_to_ndarray(c1), val.abs());
    Color {
      r: interpolated_color[0],
      g: interpolated_color[1],
      b: interpolated_color[2],
    }
  } else {
    let interpolated_color = lerp(&color_to_ndarray(black), &color_to_ndarray(c2), val.abs());
    Color {
      r: interpolated_color[0],
      g: interpolated_color[1],
      b: interpolated_color[2],
    }
  }
}

/// Creates a sparse matrix that represents the pressure relationships between cells
///
/// Rows: the number of rows in the reference grid
/// Cols: the number of columns in the reference grid
pub fn create_sprs(
  rows: usize,
  cols: usize,
) -> sprs::CsMatBase<f32, usize, Vec<usize>, Vec<usize>, Vec<f32>> {
  // we will use the triplet format for easier sparse matrix creation
  let total_cells = rows * cols;
  let mut tri = sprs::TriMat::new((total_cells, total_cells));
  for y in 0..rows {
    for x in 0..cols {
      let mut val = 0.;
      let index_total = y * cols + x;
      if x > 0 {
        val += 1.;
        tri.add_triplet(index_total, index_total - 1, -1.)
      }
      if x < cols - 1 {
        val += 1.;
        tri.add_triplet(index_total, index_total + 1, -1.)
      }
      if y > 0 {
        val += 1.;
        tri.add_triplet(index_total, index_total - cols, -1.)
      }
      if y < rows - 1 {
        val += 1.;
        tri.add_triplet(index_total, index_total + cols, -1.)
      }
      tri.add_triplet(index_total, index_total, val)
    }
  }
  tri.to_csc()
}

#[derive(Debug, Default, Clone,)]
pub struct FluidInformation {
    pub width: usize,
    pub height: usize,
    pub heat_transfer_rate: f32,
    pub bouyancy: f32,
    pub temp_max: f32,
    pub diffuse: Option<f32>,
    pub viscosity: Option<f32>,
    pub subsample: Option<u32>,
    pub confine: bool,
}
impl FluidInformation {
    /// Make a new struct with all fields set to 0 or false
    pub fn new() -> Self {
        Self {
            width: 0,
            height: 0,
            heat_transfer_rate: 0.,
            bouyancy: 0.,
            temp_max: 0.,
            diffuse: None,
            viscosity: None,
            subsample: None,
            confine: false,
        }
    }    
}
#[cfg(test)]
mod tests {
  use super::*;
  /// make a test case for the square lerp function
  #[test]
  fn test_square_lerp() {
    let point = Point { x: 0.2, y: 0.8 };
    let origin = Point { x: 0., y: 0. };
    let vals = ndarray::Array1::from_vec(vec![&3.0, &7.0, &5., &0.]);
    let result = square_lerp(point, origin, vals);
    assert_eq!(result, 5.16);
  }
  /// Tests for linear interpolation
  #[test]
  fn test_lerp() {
    let p1 = ndarray::Array1::from_vec(vec![0., 0., 0., 0.]);
    let p2 = ndarray::Array1::from_vec(vec![2., 2., 2., 2.]);
    let t = 0.5;
    let expected = ndarray::Array1::from_vec(vec![1., 1., 1., 1.]);
    let real = lerp(&p1, &p2, t);
    assert_eq!(real, expected);

    let t = 0.;
    let expected = ndarray::Array1::from_vec(vec![0., 0., 0., 0.]);
    let real = lerp(&p1, &p2, t);
    assert_eq!(real, expected);

    let t = 1.;
    let expected = ndarray::Array1::from_vec(vec![2., 2., 2., 2.]);
    let real = lerp(&p1, &p2, t);
    assert_eq!(real, expected);
  }
  // test the sparse matrix creation
  #[test]
  fn test_create_sprs() {
    let rows = 3;
    let cols = 3;
    let mat = create_sprs(rows, cols);
    assert_eq!(mat.rows(), rows * cols);
    assert_eq!(mat.cols(), rows * cols);
    for val in mat.diag_iter() {
      assert!(*val.unwrap() > 1. && *val.unwrap() < 5.);
    }
  }
}
