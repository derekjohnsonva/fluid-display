use ndarray::Array2;
use crate::utils::{self};
use sprs::CsMatBase;
// Staggered grid struct should consist of three 2d matrices
// Contents: A NxN matrix with cells that represent the contents of the fluid at that point
// X-Flow: A (N-1)xN matrix with cells that represent the x-flow at that point
// Y-Flow: A Nx(N-1) matrix with cells that represent the y-flow at that point
#[derive(Clone, Debug)]
pub struct StaggeredGrid {
    pub contents: Array2<f32>,
    pub x_flow:   Array2<f32>,
    pub y_flow:   Array2<f32>,
    // TODO: This should probably be refactored to not be in the struct
    pub sparse_rep: CsMatBase<f32, usize, Vec<usize>, Vec<usize>, Vec<f32>>, 
}
impl StaggeredGrid {
    // Constructor that takes a height and width and initializes the contents and x-flow and y-flow matrices
    pub fn new(width: usize, height: usize) -> StaggeredGrid {
        let contents =  Array2::zeros((height    , width    ));
        let x_flow =    Array2::zeros((height    , width - 1));
        let y_flow =    Array2::zeros((height - 1, width    ));
        let sparse_rep = utils::create_sprs(height, width);
        StaggeredGrid {contents, x_flow, y_flow, sparse_rep}
    }
    pub fn num_cells(&self) -> usize {
        self.contents.len()
    }
    pub fn num_rows(&self) -> &usize {
        self.contents.shape().get(0).unwrap()
    }
    pub fn num_cols(&self) -> &usize {
        self.contents.shape().get(1).unwrap()
    }
    // Make the values in the contents vector a random number between mag and -mag
    pub fn randomize_contents(&mut self, mag: f32) {
        for row in 0..self.contents.shape()[0] {
            for col in 0..self.contents.shape()[1] {
                self.contents[(row, col)] = mag * (2.0 * rand::random::<f32>() - 1.0);
            }
        }
    }
    /// Determines the force from temperature for each cell. Adds the force to the y-flow
    /// directly above and below the cell.
    pub fn add_velocity_from_temperature(&mut self, bouyancy: f32) {
        for row in 0..self.y_flow.shape()[0]  {
            for col in 0..self.y_flow.shape()[1] {
                let temp_force = self.y_flow.get((row, col)).unwrap() + bouyancy * self.get_contents(col as f32 + 0.5, row as f32 + 1.0);
                self.y_flow[(row, col)] =  temp_force;
            }
        }
    }
    // get the value of contents at a given flot point x, y
    pub fn get_contents(&self, x: f32, y: f32) -> f32 {
        // ensure that the point is within the bounds of the grid
        // the first content value is at (0.5, 0.5) on our grid and the last content value is at
        // (cols - 0.5, rows - 0.5)
        // TODO: Handle case in which row or column number is 1. Currently this breaks
        assert!(self.contents.shape()[1] > 1);
        assert!(self.contents.shape()[0] > 1);
        let max_x = self.contents.shape()[1] as f32 - 0.5;
        let max_y = self.contents.shape()[0] as f32 - 0.5;
        let x_clamped = x.clamp(0.5, max_x);
        let y_clamped = y.clamp(0.5, max_y);
        // Find the x and y vals for the content value closest to (0,0)
        let x_edge = if x_clamped < max_x {x_clamped.round()} else {x_clamped.floor()};
        let y_edge = if y_clamped < max_y {y_clamped.round()} else {y_clamped.floor()};
        let c1 = self.contents[(y_edge as usize - 1, x_edge as usize - 1)];
        let c2 = self.contents[(y_edge as usize,     x_edge as usize - 1)];
        let c3 = self.contents[(y_edge as usize - 1, x_edge as usize    )];
        let c4 = self.contents[(y_edge as usize,     x_edge as usize    )];
        let closest_points = ndarray::Array1::from(vec![&c1, &c2, &c3, &c4]);
        let point = utils::Point { x: x_clamped, y: y_clamped };
        let origin = utils::Point { x: x_edge - 0.5, y: y_edge - 0.5 };
        utils::square_lerp(point, origin, closest_points)
    }

    // Get the linearly interpolated value for x_flow at a given float x and y
    pub fn get_x_flow(&self, x: f32, y: f32) -> f32 {
        // These are the values for x and y on the nearest corner of the cell
        let x_on_corner = x.floor();
        let y_on_corner = y.round();
        // Get the x_flow values at the four verticle edges closest to the corner
        // e1 is the edge below the corner
        let e1 = self.get_x_flow_helper(x_on_corner - 1., y_on_corner - 1.);
        // e2 is the edge above the corner
        let e2 = self.get_x_flow_helper(x_on_corner - 1., y_on_corner);
        // e3 is the edge below and to the left of the corner
        let e3 = self.get_x_flow_helper(x_on_corner, y_on_corner - 1.);
        // e4 is the edge above and to the left of the corner
        let e4 = self.get_x_flow_helper(x_on_corner, y_on_corner);
        let closest_points = ndarray::Array1::from(vec![e1, e2, e3, e4]);
        let point = utils::Point { x, y };
        let origin = utils::Point { x: x_on_corner, y: y_on_corner-0.5 };
        utils::square_lerp(point, origin, closest_points)
    }
    fn get_x_flow_helper(&self, x: f32, y: f32) -> &f32 {
        if x < 0. || y < 0. {
            return &0.0
        }
        self.x_flow.get((y as usize, x as usize)).unwrap_or(&0.0)
    }
    #[allow(dead_code)]
    /// Setter for x_flow that will check if the value is within the bounds of the grid
    /// and return a result indicating if the value was set or not
    pub fn set_x_flow(&mut self, x: f32, y: f32, value: f32) -> Result<(), &str> {
        if x < 0. || x > self.x_flow.shape()[1] as f32 - 1. || y < 0. || y > self.x_flow.shape()[0] as f32 - 1. {
            return Err("Value out of bounds");
        }
        self.x_flow[(y as usize, x as usize)] = value;
        Ok(())
    }

    pub fn get_y_flow(&self, x: f32, y: f32) -> f32 {
        // These are the values for x and y on the nearest corner of the cell
        let x_on_corner = x.round();
        let y_on_corner = y.floor();
        // Get the x_flow values at the four verticle edges closest to the corner
        // e1 is the edge below the corner
        let e1 = self.get_y_flow_helper(x_on_corner - 1., y_on_corner - 1.);
        // e2 is the edge above the corner
        let e2 = self.get_y_flow_helper(x_on_corner - 1., y_on_corner);
        // e3 is the edge below and to the left of the corner
        let e3 = self.get_y_flow_helper(x_on_corner, y_on_corner - 1.);
        // e4 is the edge above and to the left of the corner
        let e4 = self.get_y_flow_helper(x_on_corner, y_on_corner);
        let closest_points = ndarray::Array1::from(vec![e1, e2, e3, e4]);
        let point = utils::Point { x, y };
        let origin = utils::Point { x: x_on_corner-0.5, y: y_on_corner};
        utils::square_lerp(point, origin, closest_points)
    }
    fn get_y_flow_helper(&self, x: f32, y: f32) -> &f32 {
        if x < 0. || y < 0. {
            return &0.0
        }
        self.y_flow.get((y as usize, x as usize)).unwrap_or(&0.0)
    }
    #[allow(dead_code)]
    /// Setter for y_flow that will check if the value is within the bounds of the grid
    /// and return a result indicating if the value was set or not
    pub fn set_y_flow(&mut self, x: usize, y: usize, value: f32) -> Result<(), &'static str> {
        if x >= self.y_flow.shape()[1] || y >= self.y_flow.shape()[0] {
            return Err("Out of bounds");
        }
        self.y_flow[(y, x)] = value;
        Ok(())
    }

    /// Get the value of the inflow in a cell
    fn divergence_helper(&self, x: usize, y: usize) -> f32 {
        // get the x-flow in the positive x direction
        let x_flow_pos = self.get_x_flow_helper(x as f32, y as f32);
        // get the x-flow in the negative x direction
        let x_flow_neg = self.get_x_flow_helper(x as f32 - 1., y as f32);
        // get the y-flow in the positive y direction
        let y_flow_pos = self.get_y_flow_helper(x as f32, y as f32);
        // get the y-flow in the negative y direction
        let y_flow_neg = self.get_y_flow_helper(x as f32, y as f32 - 1.);
        x_flow_pos - x_flow_neg + y_flow_pos - y_flow_neg
    }
    /// Find the divergence of each cell in the array and returns it as a VecViewMut. 
    /// The vector that is returned is the divergence of each cell in the grid in row major order.
    /// A positive divergence means that the flow is out of the cell.
    pub fn get_divergences(&self) -> Vec<f32> {
        let mut divergences = Vec::with_capacity(self.num_cells());
        for row in 0..self.contents.shape()[0] {
            for col in 0..self.contents.shape()[1] {
                let divergence = self.divergence_helper(col, row);
                divergences.push(divergence);
            }
        }
        divergences
    }


}
// Test the StaggeredGrid struct
#[cfg(test)]
mod tests {
    const MY_TOLERANCE: f32 = almost::F32_TOLERANCE / 5.0;

    use super::*;
    #[test]
    fn test_staggered_grid() {
        let grid = StaggeredGrid::new(10, 10);
        // Check that the contents, x-flow, and y-flow matrices have the correct size
        assert_eq!(grid.contents.shape()[0], 10);
        assert_eq!(grid.contents.shape()[1], 10);
        assert_eq!(grid.x_flow.shape()[0], 10);
        assert_eq!(grid.x_flow.shape()[1], 9);
        assert_eq!(grid.y_flow.shape()[0], 9);
        assert_eq!(grid.y_flow.shape()[1], 10);
        // Check that the contents matrix is filled with random numbers
        for i in 0..grid.contents.shape()[0] {
            for j in 0..grid.contents.shape()[1] {
                assert!(grid.contents[(i, j)] == 0.0);
            }
        }
    }
    #[test]
    fn test_staggered_grid_different_dim() {
        let rows = 5;
        let cols = 10;
        let grid = StaggeredGrid::new(cols, rows);
        // Check that the contents, x-flow, and y-flow matrices have the correct size
        assert_eq!(grid.contents.shape()[0], rows);
        assert_eq!(grid.contents.shape()[1], cols);
        assert_eq!(grid.x_flow.shape()[0], rows);
        assert_eq!(grid.x_flow.shape()[1], cols - 1);
        assert_eq!(grid.y_flow.shape()[0], rows - 1);
        assert_eq!(grid.y_flow.shape()[1], cols);
        // Check that the contents matrix is filled with random numbers
        for row in 0..grid.contents.shape()[0] {
            for col in 0..grid.contents.shape()[1] {
                assert!(grid.contents[(row, col)] == 0.0);
            }
        }
    }
    // Test the randomize_contents function
    #[test]
    fn test_randomize_contents() {
        let mut grid = StaggeredGrid::new(5, 10);
        let mag = 5.0;
        grid.randomize_contents(mag);
        // Check that the contents matrix is filled with random numbers
        grid.contents.iter().for_each(|val| {
            assert!(*val <= mag);
            assert!(*val >= -mag);
        });

        let mut grid = StaggeredGrid::new(10, 5);
        let mag = 5.0;
        grid.randomize_contents(mag);
        // Check that the contents matrix is filled with random numbers
        grid.contents.iter().for_each(|val| {
            assert!(*val <= mag);
            assert!(*val >= -mag);
        });
    }
    // Make a test for the temperature setting function
    #[test]
    fn test_add_velocity_from_temperature() {
        let mut grid = StaggeredGrid::new(2, 4);
        grid.contents[(0, 0)] = -1.0;
        grid.contents[(1, 0)] = 1.0;
        grid.contents[(2, 0)] = 1.0;
        grid.contents[(3, 0)] = 1.0;
        println!("{:?}", grid.contents);
        grid.add_velocity_from_temperature(1.0);
        println!("{}", grid.y_flow[(1,0)]);
        assert_eq!(grid.y_flow[(0, 0)], 0.0);
        assert_eq!(grid.y_flow.get((1, 0)).unwrap(), &1.0);
        assert_eq!(grid.y_flow[(2, 0)], 1.0);
    }

    #[test]
    // Test the interpolate_x_flow function
    fn test_get_x_flow() {
        let mut grid = StaggeredGrid::new(3, 3);
        grid.x_flow[(0, 0)] = 1.0;
        grid.x_flow[(1, 0)] = 1.0;
        grid.x_flow[(0, 1)] = 1.0;
        grid.x_flow[(1, 1)] = 1.0;
        let val = grid.get_x_flow(1.0, 1.0);
        assert_eq!(val, 1.0);
        let val = grid.get_x_flow(0.5, 0.5);
        assert_eq!(val, 0.5);
        let val = grid.get_x_flow(3., 3.);
        assert_eq!(val, 0.0);
        let val = grid.get_x_flow(-1., -1.);
        assert_eq!(val, 0.0);
        let val = grid.get_x_flow(-1., 1.);
        assert_eq!(val, 0.0);
    }
    #[test]
    fn test_get_x_flow_again() {
        let mut grid = StaggeredGrid::new(3, 3);
        grid.x_flow[(0, 0)] = 3.0;
        grid.x_flow[(0, 1)] = 5.0;
        grid.x_flow[(1, 0)] = 7.0;
        grid.x_flow[(1, 1)] = 0.0;
        let val = grid.get_x_flow(1.2, 1.3);
        assert!(almost::equal_with(val, 5.16, MY_TOLERANCE));

        grid.x_flow[(1, 0)] = 3.0;
        grid.x_flow[(1, 1)] = 5.0;
        grid.x_flow[(2, 0)] = 7.0;
        grid.x_flow[(2, 1)] = 0.0;
        let val = grid.get_x_flow(1.2, 2.3);
        assert!(almost::equal_with(val, 5.16, MY_TOLERANCE));
    }
    #[test]
    fn test_get_y_flow() {
        let mut grid = StaggeredGrid::new(3, 3);
        grid.y_flow[(0, 0)] = 3.0;
        grid.y_flow[(0, 1)] = 5.0;
        grid.y_flow[(1, 0)] = 7.0;
        grid.y_flow[(1, 1)] = 0.0;
        let val = grid.get_y_flow(0.7, 1.8);
        assert!(almost::equal_with(val, 5.16, MY_TOLERANCE));
    }
    #[test]
    fn test_get_contents() {
        let mut grid = StaggeredGrid::new(2, 2);
        grid.contents[(0, 0)] = 1.0;
        grid.contents[(0, 1)] = 2.0;
        grid.contents[(1, 0)] = 3.0;
        grid.contents[(1, 1)] = 4.0;
        assert_eq!(grid.get_contents(0.5, 0.5), 1.0);
        assert_eq!(grid.get_contents(0., 0.), 1.0);
        assert_eq!(grid.get_contents(0.5, 0.), 1.0);
        assert_eq!(grid.get_contents(0., 0.5), 1.0);
        assert!(almost::equal(grid.get_contents(1., 1.), 2.5));
        assert_eq!(grid.get_contents(1.5, 1.5), 4.0);
        assert_eq!(grid.get_contents(2., 1.5), 4.0);
        assert_eq!(grid.get_contents(1.5, 2.), 4.0);
        assert_eq!(grid.get_contents(2., 2.), 4.0);
        
    }
    // test find_divergences
    #[test]
    fn test_find_divergences() {
        let mut grid = StaggeredGrid::new(2, 2);
        grid.x_flow[(0, 0)] = 1.0;
        grid.x_flow[(1, 0)] = 1.0;
        grid.y_flow[(0, 0)] = 1.0;
        grid.y_flow[(0, 1)] = 1.0;
        let divergences = grid.get_divergences();
        assert_eq!(divergences[0], 2.0);
        assert_eq!(divergences[1], 0.0);
        assert_eq!(divergences[2], 0.0);
        assert_eq!(divergences[3], -2.0);
    }
    #[test]
    fn test_find_divergences2() {
        let mut grid = StaggeredGrid::new(2, 2);
        //    | 3 | 4 |
        //    | 1 | 2 |
        // (0,0)
        grid.x_flow[(0, 0)] = 1.0;
        grid.x_flow[(1, 0)] = 0.0;
        grid.y_flow[(0, 0)] = 1.0;
        grid.y_flow[(0, 1)] = 4.0;
        let divergences = grid.get_divergences();
        assert_eq!(divergences[0], 2.0);
        assert_eq!(divergences[1], 3.0);
        assert_eq!(divergences[2], -1.0);
        assert_eq!(divergences[3], -4.0);
    }
}