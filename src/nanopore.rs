//! # Flowcell and Condition utilities.
//!
//! This module provides utility functions for working with flowcells. It includes functions for retrieving the coordinates of channels on a flowcell (`get_coords`), generating a flowcell layout (`get_flowcell_array`), and dividing the flowcell into sections (`generate_flowcell`).
//!
//! ## Examples
//!
//! ```rust, ignore
//! use flowcell_utils::{get_coords, generate_flowcell, get_flowcell_array};
//!
//! // Retrieve the coordinates of a channel on a flowcell
//! let coords = get_coords(5, 126);
//! println!("Channel 5 coordinates: {:?}", coords);
//!
//! // Generate the layout of a flowcell
//! let flowcell_layout = get_flowcell_array(512);
//! println!("Flowcell layout: {:?}", flowcell_layout);
//!
//! // Divide the flowcell into sections
//! let divided_flowcell = generate_flowcell(512, 4, 0, false);
//! println!("Divided flowcell: {:?}", divided_flowcell);
//! ```
//!
use crate::channels::{FLONGLE_CHANNELS, MINION_CHANNELS};
use ndarray::{s, Array, Array2, Axis};
/// Returns the coordinates (column, row) of a channel on a flowcell.
///
/// # Arguments
///
/// * `channel` - The channel to retrieve the coordinates for.
/// * `flowcell_size` - The total number of channels on the flowcell.
///
/// # Returns
///
/// * `Ok((column, row))` - The column and row of the channel number in the flowcell.
/// * `Err(error_message)` - An error message indicating the reason for failure.
///
/// # Errors
///
/// This function may return an error in the following cases:
///
/// * If the `channel` is greater than the `flowcell_size`.
/// * If the `flowcell_size` is not recognized (not equal to 3000, 126, or 512).
/// * If the `channel` is not found in the predefined coordinate dictionaries for Flongle or MinION channels.
///
/// # Examples
///
/// ```
/// use readfish_tools::nanopore::get_coords;
///
/// let result = get_coords(5, 126);
/// assert_eq!(result, Ok((4, 9)));
///
/// let result = get_coords(300, 3000);
/// assert_eq!(result, Ok((19, 4)));
///
/// let result = get_coords(0, 512);
/// assert_eq!(result, Err("channel cannot be less than one or above flowcell_size".to_string()));
///
/// let result = get_coords(1, 513);
/// assert_eq!(result, Err("flowcell_size is not recognized".to_string()));
/// ```
pub fn get_coords(channel: usize, flowcell_size: usize) -> Result<(usize, usize), String> {
    if (channel > flowcell_size) | (channel == 0) {
        return Err("channel cannot be less than one or above flowcell_size".to_string());
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

/// Returns an `Array2` representing the layout of a flowcell.
///
/// The flowcell layout is generated based on the provided `flowcell_size`. Each channel is assigned coordinates (column, row)
/// using the `get_coords` function. The maximum row and column values are determined from the channel coordinates, and an `Array2`
/// is initialized with these dimensions. The flowcell layout is then mimicked in the array by adding the channel numbers in the
/// respective coordinates. The resulting array rows are reversed to obtain the correct orientation.
///
/// # Arguments
///
/// * `flowcell_size` - The total number of channels on the flowcell.
///
/// # Returns
///
/// An `Array2` representing the layout of the flowcell.
///
/// # Panics
///
/// This function may panic if the `get_coords` function returns an error.
///
/// # Examples
///
/// ```rust,ignore
/// use crate::get_flowcell_array;
/// use ndarray::array;
///
/// let result = get_flowcell_array(512);
/// // [[121,113,...], [122, 114,...],...]
///
///
/// ```
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

/// Generates a flowcell divided into sections based on the provided parameters.
///
/// If `odd_even` is `true`, the function returns a vector of two vectors where the first vector contains odd channels
/// and the second vector contains even channels, from 1 to `flowcell_size`.
///
/// If `odd_even` is `false`, the function generates a 2D `Array2` representing the layout of the flowcell using the
/// `get_flowcell_array` function. The flowcell is then divided into `split` sections along the specified `axis` (0 for
/// rows, 1 for columns). The number of sections must evenly divide the length of the target `axis` dimension. The
/// resulting divided flowcell is returned as a vector of vectors.
///
///
///
/// # Arguments
///
/// * `flowcell_size` - The total number of channels on the flowcell.
/// * `split` - The number of sections to split the flowcell into.
/// * `axis` - The axis along which to split the flowcell (0 for rows, 1 for columns).
/// * `odd_even` - Specifies whether to return the flowcell divided into odd and even channels.
///
/// # Panics
///
/// This function may panic in the following cases:
///
/// * If `split` is 0, indicating an invalid value for the number of sections.
/// * If the target axis dimension cannot be evenly divided by `split`, resulting in an uneven split.
///
/// # Examples
///
/// ```rust,ignore
/// use crate::{generate_flowcell, get_flowcell_array};
/// use ndarray::array;
///
/// generate_flowcell(512, 4, 1, false);
/// // splits the flowcell into vertical quarters
/// // ########........................
/// // ########........................
/// // ########........................
/// // ########........................
///
/// // moving across the flowcell quarters, where # is the a channel in this split.
/// generate_flowcell(512, 4, 1, false);
/// // [
/// //    [121, 113, 105, 97, 185, 177, 169, 161, ... 129],
/// //    [249, 241, 233, 225, 313, 305, 297, 289, ... 257],
/// //    [377, 369, 361, 353, 441, 433, 425, 417, ... 385],
/// //    [505, 497, 489, 481, 57, 49, 41, 33, 506, ... 1]
/// // ]
///
/// generate_flowcell(512, 2, 1, true);
/// // [
/// //    [1, 3, 5, 7, 9, 11, 13,...
/// //    [2, 4, 6, 8, 10, 12, 14,...
/// // ]
pub fn generate_flowcell(
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

    split_flowcell
}

/// Formats a given number of bases into a human-readable string with appropriate units (Kb, Mb, Gb, etc.).
///
/// # Arguments
///
/// * `number` - The number of bases to be formatted.
///
/// # Returns
///
/// A string representing the formatted number of bases with the appropriate unit.
///
/// # Examples
///
/// ```rust,ignore
/// assert_eq!(format_bases(1_000), "1.00 Kb");
/// assert_eq!(format_bases(1_000_000), "1.00 Mb");
/// assert_eq!(format_bases(1_630_000), "1.63 Mb");
/// assert_eq!(format_bases(1_000_000_000), "1.00 Gb");
/// ```
pub fn format_bases(number: usize) -> String {
    let number = number as f64;
    let units = ["", "K", "M", "G", "T", "P", "E", "Z", "Y"];
    let base = 1000.0;

    if number.abs() < base {
        return format!("{} b", number);
    }

    let exponent = (number.abs().log(base)).floor() as i32;
    let unit_idx = if exponent >= 0 {
        std::cmp::min(exponent as usize, units.len() - 1)
    } else {
        std::cmp::max(exponent as usize + units.len(), 0)
    };

    let formatted_number = number / base.powi(exponent);
    format!("{:.2} {}b", formatted_number, units[unit_idx])
}

/// Calculate the running mean incrementally.
///
/// Given a mutable reference to `mean`, `count`, and `value`, this function calculates the running
/// mean by updating the variables in place.
///
/// # Arguments
///
/// * `mean`: A mutable reference to the current running mean.
/// * `count`: A mutable reference to the count of elements seen so far.
/// * `value`: A mutable reference to the new value to be included in the running mean calculation.
///
/// # Example
///
/// ```
/// use readfish_tools::nanopore::running_mean;
/// let mut mean = 0;
/// let mut count = 0;
/// let mut value = 5;
/// running_mean(&mut mean, &mut count, &mut value);
/// assert_eq!(mean, 5);
/// assert_eq!(count, 1);
/// ```
pub fn running_mean(mean: &mut isize, count: &mut isize, value: &mut isize) {
    *count += 1;
    *mean += (*value - *mean) / *count; // Update the running mean incrementally
}

// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_running_mean() {
        let mut mean = 0;
        let mut count = 0;
        let mut value = 5;

        running_mean(&mut mean, &mut count, &mut value);
        assert_eq!(mean, 5);
        assert_eq!(count, 1);

        // Add more values and update the running mean
        value = 8;
        running_mean(&mut mean, &mut count, &mut value);
        assert_eq!(mean, 6); // (5 + 8) / 2 = 6
        assert_eq!(count, 2);

        value = 12;
        running_mean(&mut mean, &mut count, &mut value);
        assert_eq!(mean, 8); // (5 + 8 + 12) / 3 = 8
        assert_eq!(count, 3);

        value = 4;
        running_mean(&mut mean, &mut count, &mut value);
        assert_eq!(mean, 7); // (5 + 8 + 12 + 4) / 4 = 7
        assert_eq!(count, 4);
    }

    #[test]
    fn test_format_bases() {
        assert_eq!(format_bases(1_000), "1.00 Kb");
        assert_eq!(format_bases(1_000_000), "1.00 Mb");
        assert_eq!(format_bases(1_630_000), "1.63 Mb");
        assert_eq!(format_bases(1_000_000_000), "1.00 Gb");
    }
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
        assert_eq!(fa.get((0, 0)).unwrap(), &121_usize)
    }

    #[test]
    #[should_panic]
    fn test_get_flowcell_array_panic() {
        let fa = get_flowcell_array(513);
        assert_eq!(fa.get((0, 0)).unwrap(), &121_usize)
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
