//!Hello
//!
use crate::channels::{FLONGLE_CHANNELS, MINION_CHANNELS};

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

// Tests
#[cfg(test)]
mod tests {
    use super::*;

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
