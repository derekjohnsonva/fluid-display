use ndarray::Array2;
use crate::utils::FluidInformation;
use crate::{diffusion, staggered_grid::StaggeredGrid};

/// Takes in a reference to a ndarray, sets each value in the bottom row to
/// h + (1 - h) * current_value
pub fn update_bottom_row(arr: &mut Array2<f32>, h: f32) {
    for col in 0..arr.shape()[1] {
        arr[(0, col)] = h + (1.0 - h) * arr[(0, col)];
    }
}
/// Takes in a reference to a ndarray, sets each value in the top row to
/// -h + (1 - h) * current_value
pub fn update_top_row(arr: &mut Array2<f32>, h: f32) {
    let top_row = arr.shape()[0] - 1;
    for col in 0..arr.shape()[1] {
        arr[(top_row, col)] = -h + (1.0 - h) * arr[(top_row, col)];
    }
}
/// Takes in a source grid and a destination grid. Updates the destination grid
/// to have contents back advected from the source grid.
///
/// Info on back advection: [source](https://05f43270-a325-4aeb-a1d8-39e42f0e9fcc.filesusr.com/ugd/cf1fd6_ad9e7e5ca0c24e4aae415c84f8217e10.pdf)
pub fn advect_content(dest_grid: &mut StaggeredGrid, source_grid: &StaggeredGrid, dt: f32) {
    (0..*dest_grid.num_rows()).for_each(|row| {
        (0..*dest_grid.num_cols()).for_each(|col| {
            // find the x and y velocity at a the center of the current cell
            let x_vel = dest_grid.get_x_flow(col as f32 + 0.5, row as f32 + 0.5);
            let y_vel = dest_grid.get_y_flow(col as f32 + 0.5, row as f32 + 0.5);
            let dx = x_vel * dt;
            let new_x = col as f32 + 0.5 - dx;
            let dy = y_vel * dt;
            let new_y = row as f32 + 0.5 - dy;
            dest_grid.contents[(row, col)] = source_grid.get_contents(new_x, new_y);
        });
    });
}
/// Takes in a source grid and a destination grid. Updates the destination grid
/// to have y-flow back advected from the source grid.
pub fn advect_x_flow(dest_grid: &mut StaggeredGrid, source_grid: &StaggeredGrid, dt: f32) {
    for row in 0..dest_grid.x_flow.shape()[0] {
        for col in 0..dest_grid.x_flow.shape()[1] {
            // find the x and y velocity at a the center of the current cell
            let x_actual = col as f32 + 1.0;
            let y_actual = row as f32 + 0.5;
            let x_vel = source_grid.get_x_flow(x_actual, y_actual);
            let y_vel = source_grid.get_y_flow(x_actual, y_actual);
            let dx = x_vel * dt;
            let new_x = x_actual - dx;
            let dy = y_vel * dt;
            let new_y = y_actual - dy;
            dest_grid.x_flow[(row, col)] = source_grid.get_x_flow(new_x, new_y);
        }
    }
}
/// Takes in a source grid and a destination grid. Updates the destination grid
/// to have y-flow back advected from the source grid.
pub fn advect_y_flow(dest_grid: &mut StaggeredGrid, source_grid: &StaggeredGrid, dt: f32) {
    for row in 0..dest_grid.y_flow.shape()[0] {
        for col in 0..dest_grid.y_flow.shape()[1] {
            // find the x and y velocity at a the center of the current cell
            let x_actual = col as f32 + 0.5;
            let y_actual = row as f32 + 1.0;
            let x_vel = source_grid.get_x_flow(x_actual, y_actual);
            let y_vel = source_grid.get_y_flow(x_actual, y_actual);
            let dx = x_vel * dt;
            let new_x = x_actual - dx;
            let dy = y_vel * dt;
            let new_y = y_actual - dy;
            dest_grid.y_flow[(row, col)] = source_grid.get_y_flow(new_x, new_y);
        }
    }
}
/// Performs one cycle of simulation on the current grid.
/// A cycle consists of the following steps:
/// 1. Add temperature values to the top and bottom rows of the grid
/// 2. Update the velocities based on the current temperature values
/// 3. Diffuse the content, x-flow, and y-flow of the grid
/// 4. Advect the x-flow and y-flow
/// 5. Advect the content of the grid
pub fn simulate_step(grid: &mut StaggeredGrid, dt: f32, fluid_info: &FluidInformation) {
    let h = fluid_info.heat_transfer_rate;
    let bouyancy = fluid_info.bouyancy;
    add_temp_velocity(grid, h, bouyancy);
    let mut past_grid = grid.clone();
    // get the sum of the array
    diffuse(grid, &mut past_grid, fluid_info);
    // get the sum of the past_grid content
    // This has to come after diffusion because diffusion is
    // supposed to change the flow and content totals
    // TODO: Determine a better way to do this
    let mut x_flow_sum_before = 0.0;
    let mut y_flow_sum_before = 0.0;
    if fluid_info.confine {
        grid.x_flow.iter().for_each(|f| x_flow_sum_before += f.powi(2));
        grid.y_flow.iter().for_each(|f| y_flow_sum_before += f.powi(2));
    }
    vel_step(grid, &mut past_grid, dt);
    dens_step(grid, &mut past_grid, dt);
    if fluid_info.confine {
        let mut x_flow_sum_after = 0.0;
        let mut y_flow_sum_after = 0.0;
        for i in 0..grid.x_flow.shape()[0] {
            for j in 0..grid.x_flow.shape()[1] {
                x_flow_sum_after += grid.x_flow.get((i, j)).unwrap().powi(2);
            }
        }
        for i in 0..grid.y_flow.shape()[0] {
            for j in 0..grid.y_flow.shape()[1] {
                y_flow_sum_after += grid.y_flow.get((i, j)).unwrap().powi(2);
            }
        }
        let x_flow_mult = (x_flow_sum_before / x_flow_sum_after).sqrt();
        let y_flow_mult = (y_flow_sum_before / y_flow_sum_after).sqrt();
        for i in 0..grid.x_flow.shape()[0] {
            for j in 0..grid.x_flow.shape()[1] {
                let new_val = grid.x_flow.get((i, j)).unwrap() * x_flow_mult;
                grid.x_flow[(i, j)] =  new_val;
            }
        }
        for i in 0..grid.y_flow.shape()[0] {
            for j in 0..grid.y_flow.shape()[1] {
                let new_val = grid.y_flow.get((i, j)).unwrap() * y_flow_mult;
                grid.y_flow[(i, j)] = new_val;
            }
        }
    }
}
pub fn add_temp_velocity(grid: &mut StaggeredGrid, h: f32, bouyancy: f32) {
    update_bottom_row(&mut grid.contents, h); // add heat to the bottom row
    update_top_row(&mut grid.contents, h); // add cooling the the top row
    grid.add_velocity_from_temperature(bouyancy);
    project(grid);
}

fn diffuse(
    source_grid: &mut StaggeredGrid,
    dest_grid: &mut StaggeredGrid,
    fluid_info: &FluidInformation,
) {
    if let Some(diffuse) = fluid_info.diffuse {
        //TODO: make filter size a global constant or a param
        diffusion::filter_array(
            &mut source_grid.contents,
            &mut dest_grid.contents,
            diffuse,
            5,
            Some(10e-6),
        );
    }

    if let Some(diffuse) = fluid_info.viscosity {
        println!("Diffusion!!!");
        diffusion::filter_array(
            &mut source_grid.x_flow,
            &mut dest_grid.x_flow,
            diffuse,
            7,
            Some(10e-6),
        );
        diffusion::filter_array(
            &mut source_grid.y_flow,
            &mut dest_grid.y_flow,
            diffuse,
            7,
            Some(10e-6),
        );
    }

    project(dest_grid);
    // TODO: Not sure if this is needed
    *source_grid = dest_grid.clone();
}

pub fn dens_step(cur_grid: &mut StaggeredGrid, past_grid: &mut StaggeredGrid, dt: f32) {
    advect_content(cur_grid, past_grid, dt);
}

pub fn vel_step(cur_grid: &mut StaggeredGrid, past_grid: &mut StaggeredGrid, dt: f32) {
    advect_x_flow(cur_grid, past_grid, dt);
    advect_y_flow(cur_grid, past_grid, dt);
    project(cur_grid);
}

fn project(cur_grid: &mut StaggeredGrid) {
    let (gradient_x, gradient_y) = make_pressure_gradient(cur_grid);
    cur_grid.x_flow = &cur_grid.x_flow - gradient_x;
    cur_grid.y_flow = &cur_grid.y_flow - gradient_y;
}

pub fn make_pressure_gradient(grid: &StaggeredGrid) -> (Array2<f32>, Array2<f32>) {
    let divergences = OwnedVec::from(grid.get_divergences());
    let mut pressure = OwnedVec::from_elem(divergences.shape()[0], 0.0);
    let (_iter, _error) = iterative_solver(
        grid.sparse_rep.view(),
        pressure.view_mut(),
        divergences.view(),
        20,
        10e-10,
    );
    // now that we have found pressure, we can find the velocities for the
    // incompressible field
    let mut x_flow_gradient =
        Array2::zeros(grid.x_flow.raw_dim());
    let mut y_flow_gradient =
        Array2::zeros(grid.y_flow.raw_dim());
    for row in 0..*grid.num_rows() {
        for col in 0..*grid.num_cols() {
            // iterate through each cell and find the x-vel to the right and the y-vel down
            if col < grid.num_cols() - 1 {
                let pressure_index = row * grid.num_cols() + col;
                let mut pressure_delta = pressure[pressure_index] - pressure[pressure_index + 1];
                if pressure_delta.is_nan() {
                    pressure_delta = 0.0
                }
                x_flow_gradient[(row, col)] = pressure_delta;
            }
            if row < grid.num_rows() - 1 {
                let pressure_index = row * grid.num_cols() + col;
                let mut pressure_delta =
                    pressure[pressure_index] - pressure[pressure_index + grid.num_cols()];
                if pressure_delta.is_nan() {
                    pressure_delta = 0.0
                }
                y_flow_gradient[(row, col)] = pressure_delta;
            }
        }
    }
    (x_flow_gradient, y_flow_gradient)
}

// *************************Sparse Matrix Solving************************************** //
type VecView<'a, T> = ndarray::ArrayView<'a, T, ndarray::Ix1>;
type VecViewMut<'a, T> = ndarray::ArrayViewMut<'a, T, ndarray::Ix1>;
type OwnedVec<T> = ndarray::Array<T, ndarray::Ix1>;

fn iterative_solver(
    mat: sprs::CsMatView<f32>,
    mut x: VecViewMut<f32>,
    rhs: VecView<f32>,
    max_iter: usize,
    error_goal: f32,
) -> (usize, f32) {
    assert!(mat.rows() == mat.cols());
    assert!(mat.rows() == x.shape()[0]);
    let mut r = &rhs - (&mat * &x); // the residual vector; stores the error of the current estimate
    let mut x_movement = r.clone(); // The direction we'll move x in to improve it
    let mut error = r.dot(&r); // the squared magnitude of r
    for it in 0..max_iter {
        let x_move = x_movement; // temporary value
        let alpah_p = &mat * &x_move;
        let alpha = error / x_move.dot(&alpah_p); // how far to move in the x_movement direction
        x += &(alpha * &x_move);
        r = r - alpha * &alpah_p; // update the error
        let new_error = r.dot(&r);
        let beta = new_error / error; // used to ensure each direction is conjugate with all previous directions
        error = new_error; // update the error squared magnitude
        x_movement = x_move * beta; // update the direction for the next step
        x_movement += &r; // update the direction for the next step
                                      // error corresponds to the state before iteration, but
                                      // that shouldn't be a problem
        if error < error_goal {
            return (it, error);
        }
    }
    (max_iter - 1, error)
}
#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;
    const MY_TOLERANCE: f32 = almost::F32_TOLERANCE;
    #[test]
    fn test_update_bottom_row() {
        let mut arr = Array2::from_elem((3, 3), 0.0);
        arr[(0, 0)] = 1.0;
        arr[(0, 1)] = 2.0;
        arr[(0, 2)] = 3.0;
        update_bottom_row(&mut arr, 0.5);
        assert_eq!(arr[(0, 0)], 1.0);
        assert_eq!(arr[(0, 1)], 1.5);
        assert_eq!(arr[(0, 2)], 2.);
        update_bottom_row(&mut arr, 0.5);
        assert_eq!(arr[(0, 0)], 1.0);
        assert_eq!(arr[(0, 1)], 1.25);
        assert_eq!(arr[(0, 2)], 1.5);
    }

    #[test]
    fn test_contect_advection() {
        let mut grid = StaggeredGrid::new(2, 2);
        grid.contents[(0, 0)] = 1.0;
        grid.contents[(0, 1)] = 2.0;
        grid.contents[(1, 0)] = 3.0;
        grid.contents[(1, 1)] = 4.0;
        grid.x_flow[(0, 0)] = -1.;
        grid.x_flow[(1, 0)] = -1.;
        grid.y_flow[(0, 0)] = -1.;
        grid.y_flow[(0, 1)] = -1.;
        let mut grid_past = grid.clone();
        let dt = 2.;
        advect_content(&mut grid, &mut grid_past, dt);
        println!("{}", grid.contents[(0, 0)]);
        assert!(almost::equal_with(grid.contents[(0, 0)], 4.0, MY_TOLERANCE));
    }

    #[test]
    fn test_gauss_seidel() {
        let trip = TriMat::from_triplets(
            (3, 3),
            vec![0, 0, 0, 1, 1, 1, 2, 2, 2],
            vec![0, 1, 2, 0, 1, 2, 0, 1, 2],
            vec![4., -1., -1., -2., 6., 1., -1., 1., 7.],
        );
        let mat = trip.to_csc();
        let mut rhs = OwnedVec::from_elem(3, 1.0);
        rhs[[0]] = 3.0;
        rhs[[1]] = 9.0;
        rhs[[2]] = -6.0;
        let mut x = OwnedVec::from_elem(3, 0.0);
        let goal_error = 1e-10;
        let max_iter = 100;
        let (it, error) =
            iterative_solver(mat.view(), x.view_mut(), rhs.view(), max_iter, goal_error);
        assert!(it < max_iter);
        assert!(error <= goal_error);
        assert!(almost::equal_with(x[[0]], 1., MY_TOLERANCE)); // should be 1
        assert!(almost::equal_with(x[[1]], 2., MY_TOLERANCE)); // should be 2
        assert!(almost::equal_with(x[[2]], -1., MY_TOLERANCE)); // should be -1
    }
    #[test]
    fn test_make_pressure_gradient() {
        let grid = StaggeredGrid::new(2, 2);
        let (x_flow_incompressible, y_flow_incompressible) = make_pressure_gradient(&grid);
        x_flow_incompressible
            .iter()
            .for_each(|val| {
                assert_eq!(*val, 0.0);
            });
        y_flow_incompressible
            .iter()
            .for_each(|val| {
                assert_eq!(*val, 0.0);
            });
    }
    #[test]
    fn test_make_pressure_gradient2() {
        let mut grid = StaggeredGrid::new(2, 2);
 
        grid.x_flow[(0, 0)] =  1.;
        grid.x_flow[(1, 0)] =  1.;
        grid.y_flow[(0, 0)] =  1.;
        grid.y_flow[(0, 1)] =  1.;
        let (x_flow_incompressible, y_flow_incompressible) = make_pressure_gradient(&grid);
        assert_eq!(x_flow_incompressible[(0, 0)], 1.0);
        assert_eq!(x_flow_incompressible[(1, 0)], 1.0);
        assert_eq!(y_flow_incompressible[(0, 0)], 1.0);
        assert_eq!(y_flow_incompressible[(0, 1)], 1.0);

        grid.x_flow[(0, 0)] =  1.;
        grid.x_flow[(1, 0)] =  0.;
        grid.y_flow[(0, 0)] =  1.;
        grid.y_flow[(0, 1)] =  4.;
        let (x_flow_incompressible, y_flow_incompressible) = make_pressure_gradient(&grid);
        assert_eq!(x_flow_incompressible[(0, 0)], 0.0);
        assert_eq!(x_flow_incompressible[(1, 0)], 1.0);
        assert_eq!(y_flow_incompressible[(0, 0)], 2.0);
        assert_eq!(y_flow_incompressible[(0, 1)], 3.0);
    }
}
