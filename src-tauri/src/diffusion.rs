use ndarray::Array2;
fn filter_element(n: usize, mut t: f32, cutoff: Option<f32>) -> f32 {
    let mut answer = 0.0;
    let mut term = 1.0;
    let scale = (-t).exp();
    t /= 2.;
    for i in 1..n+1 {term *= t/(i as f32)}
    t *= t;
    let mut m = 0.0;
    let cut = cutoff.unwrap_or(0.0);
    while term.abs() > cut {
        answer += term;
        m += 1.0;
        term *= t / (m *(m + n as f32));
    }
    scale * answer
}

fn make_filter(size: usize, t: f32, cutoff: Option<f32>) -> Vec<f32> {
    assert!(size % 2 == 1); // filter should be an odd size
    let center_val = filter_element(0, t, cutoff);
    let mut outer_vals = Vec::new();
    let outer_indexes = size / 2;
    for i in 1..outer_indexes+1 {
        outer_vals.push(filter_element(i, t, cutoff))
    }
    let mut front_vec = outer_vals.clone();
    front_vec.reverse();
    front_vec.push(center_val);
    front_vec.append(&mut outer_vals);
    assert_eq!(front_vec.len(), size);
    front_vec
}

fn get_filter_row(arr: &mut Array2<f32>, row: usize, col: usize, size: usize) -> Vec<&f32> {
    assert!(size % 2 == 1);
    let mut row_vec = Vec::with_capacity(size);
    let col_indexes: Vec<usize> = ((col as i32 - (size as i32)/2)..(col as i32 + (size as i32)/2) + 1).map(|x| x.clamp(0, arr.shape()[1] as i32-1) as usize).collect();   
    for col_index in col_indexes {
        let val = arr.get((row, col_index)).unwrap_or(&0.0);
        row_vec.push(val);
    }
    assert_eq!(row_vec.len(), size);
    row_vec

}
fn get_filter_col(arr: &mut Array2<f32>, row: usize, col: usize, size: usize) -> Vec<&f32> {
    assert!(size % 2 == 1);
    let mut col_vec = Vec::with_capacity(size);
    let row_indexes: Vec<usize> = ((row as i32 - (size as i32)/2)..(row as i32 + (size as i32)/2) + 1).map(|x| x.clamp(0, arr.shape()[0] as i32-1) as usize).collect();   
    for row_index in row_indexes {
        let val = arr.get((row_index, col)).unwrap_or(&0.0);
        col_vec.push(val);
    }
    assert_eq!(col_vec.len(), size);
    col_vec
}

pub fn filter_array(source_arr: &mut Array2<f32>, dest_array: &mut Array2<f32>, t: f32, filter_size: usize, cutoff: Option<f32>) {
    // TODO: This should only be made once
    let filter = make_filter(filter_size, t, cutoff);
    for row in 0..source_arr.shape()[0] {
        for col in 0..source_arr.shape()[1] {
            let row_vec = get_filter_row(source_arr, row, col, filter_size);
            // multiply the values of row vec and filter together and get the sum
            let mut sum = 0.0;
            for i in 0..filter.len() {
                sum += row_vec[i] * filter[i];
            }
            // assert!(sum.abs() <= 1.0);
            dest_array[(row, col)] =  sum.clamp(-1., 1.);
        }
    }
    for row in 0..source_arr.shape()[0] {
        for col in 0..source_arr.shape()[1] {
            let col_vec = get_filter_col(source_arr, row, col, filter_size);
            // multiply the values of row vec and filter together and get the sum
            let mut sum = 0.0;
            for i in 0..filter.len() {
                sum += col_vec[i] * filter[i];
            }
            // assert!(sum.abs() <= 1.0);
            dest_array[(row, col)] =  sum.clamp(-1., 1.0);
        }
    }
}


// test module
#[cfg(test)]
mod tests {
    use super::*;
    const MY_TOLERANCE: f32 = almost::F32_TOLERANCE;
    #[test]
    fn test_make_filter() {
        let filter = make_filter(15, 1.5, None);
        let sum = filter.iter().fold(0.0, |acc, x| acc + x);
        assert!(almost::equal_with(sum, 1.0, MY_TOLERANCE));

        let filter = make_filter(15, 1.5, Some(10e-6));
        let sum = filter.iter().fold(0.0, |acc, x| acc + x);
        assert!(almost::equal_with(sum, 1.0, MY_TOLERANCE));
        let filter = make_filter(9, 1.2, Some(10e-6));
        let sum = filter.iter().fold(0.0, |acc, x| acc + x);
        assert!(sum < 1.0);
    }
    #[test]
    fn test_filter_array() {
        let mut arr = Array2::from_elem((20, 20), 1.0);
        let mut dest_arr = Array2::from_elem((20, 20), 0.0);
        filter_array(&mut arr, &mut dest_arr, 1.5, 15, None);
        for val in dest_arr.iter() {
            assert!(almost::equal_with(*val, 1.0, MY_TOLERANCE));
        }
        let mut arr = Array2::from_elem((20, 20), 0.0);
        let mut dest_arr = Array2::from_elem((20, 20), 1.0);
        filter_array(&mut arr, &mut dest_arr, 1.5, 5, None);
        for val in dest_arr.iter() {
            assert!(almost::equal_with(*val, 0.0, MY_TOLERANCE));
        }

    }
}