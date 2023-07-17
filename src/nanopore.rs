//!Hello
//!
use crate::channels::{FLONGLE_CHANNELS, MINION_CHANNELS};
use ndarray::{s, Array, Array2, ArrayView2, Axis};
/// Get the cartesian coordinates for a given channel on the flowcell
pub fn get_coords(channel: usize, flowcell_size: usize) -> Result<(usize, usize), String> {
    if channel > flowcell_size {
        return Err("channel cannot be below 0 or above flowcell_size".to_string());
    }

    if flowcell_size == 3000 {
        // find which block of 12 we are in
        let block = (channel - 1) / 250;
        let remainder = (channel - 1) % 250;
        let row = remainder / 10;
        let column = remainder % 10 + block * 10;
        Ok((column, row))
    } else if flowcell_size == 126 {
        match FLONGLE_CHANNELS.get(&channel) {
            Some(coordinates) => Ok(*coordinates),
            None => Err("channel not found in FLONGLE_CHANNELS".to_string()),
        }
    } else if flowcell_size == 512 {
        match MINION_CHANNELS.get(&channel) {
            Some(coordinates) => Ok(*coordinates),
            None => Err("channel not found in MINION_CHANNELS".to_string()),
        }
    } else {
        Err("flowcell_size is not recognized".to_string())
    }
}

/// Get a 2d array representation of an array, with the channel numbers placed correctly
fn get_flowcell_array(flowcell_size: usize) -> Array2<usize> {
    // Make a vector of tuples of (column, row, channel)
    let coords: Vec<(usize, usize, usize)> = (1..=flowcell_size)
        .map(|x| {
            let (col, row) = get_coords(x, flowcell_size).unwrap();
            (col, row, x)
        })
        .collect();

    // Determine the maximum row and column from the coords vector
    let max_row = coords.iter().map(|&(_, row, _)| row).max().unwrap();
    let max_col = coords.iter().map(|&(col, _, _)| col).max().unwrap();

    // Initialize an Array2 using the max row and column values
    let mut flowcell_layout = Array::zeros((max_row + 1, max_col + 1));

    // Mimic flowcell layout in an array
    for &(col, row, chan) in &coords {
        flowcell_layout[[row, col]] += chan;
    }

    // return the reversed array, to get the right orientation
    flowcell_layout.slice(s![..;-1,..]).to_owned()
}

fn generate_flowcell(
    flowcell_size: usize,
    split: usize,
    axis: usize,
    odd_even: bool,
) -> Vec<Vec<usize>> {
    if odd_even {
        return vec![
            (1..=flowcell_size).step_by(2).collect(),
            (2..=flowcell_size).step_by(2).collect(),
        ];
    }

    let arr: Array2<usize> = get_flowcell_array(flowcell_size);

    if split == 0 {
        panic!("split must be a positive integer");
    }

    let (dim1, dim2) = arr.dim();
    let target_dim = if axis == 0 { dim1 } else { dim2 };

    if target_dim % split != 0 {
        panic!("The flowcell cannot be split evenly");
    }
    let axis_ = Axis(axis);
    let split_flowcell = arr
        .axis_chunks_iter(axis_, arr.len_of(axis_) / split)
        .map(|x| x.iter().cloned().collect())
        .collect::<Vec<Vec<usize>>>();

    // let chunk_size = target_dim / split;

    // let mut split_arr = if axis == 0 {
    //     arr.axis_chunks_iter(Axis(0), chunk_size)
    // } else {
    //     arr.axis_chunks_iter(Axis(1), chunk_size)
    // };

    // let mut result = Vec::new();
    // for _ in 0..split {
    //     let chunk: ArrayView2<usize> = split_arr.next().unwrap();
    //     result.push(chunk.iter().cloned().collect());
    // }

    split_flowcell
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_flowcell() {
        let x = generate_flowcell(512, 2, 1, false);
        assert_eq!(x.len(), 2);
        assert_eq!(x[0][0], 121_usize);
        assert_eq!(x[1][0], 377_usize)
    }

    #[test]
    fn test_generate_flowcell_odd_even() {
        let x = generate_flowcell(512, 0, 0, true);
        assert_eq!(x.len(), 2);
        assert_eq!(x[0][0], 1);
        assert_eq!(x[1][0], 2)
    }

    #[test]
    fn test_get_flowcell_array() {
        let fa = get_flowcell_array(512);
        assert_eq!(fa.get((0, 0)).unwrap(), &89_usize)
    }

    #[test]
    #[should_panic]
    fn test_get_flowcell_array_panic() {
        let fa = get_flowcell_array(513);
        assert_eq!(fa.get((0, 0)).unwrap(), &89_usize)
    }

    #[test]
    fn test_get_coords() {
        assert_eq!(get_coords(2, 512).unwrap(), (31_usize, 1_usize));
        assert_eq!(get_coords(2, 126).unwrap(), (1_usize, 9_usize));
        assert_eq!(get_coords(2, 3000).unwrap(), (1_usize, 0_usize));
    }

    #[test]
    #[should_panic]
    fn test_get_coords_panics() {
        // Code that is expected to panic
        get_coords(10000, 10).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_get_coords_panics_size() {
        // Code that is expected to panic
        get_coords(10, 127).unwrap();
    }
}
