//! Manipulate a sequencing summary file.
//! Currently we cannot demultiplex the sequencing summary file, so we need to manipulate the sequencing summary file
//! to match the demultiplexed PAF file on Read ID. This is the purpose of the `get_record` function.
//!
//! A buffer is used to store the sequencing summary records. The buffer is a linked hash map,
//!  with read ID as the key and tuples containing `SeqSumInfo` variants as the values.
//! Currently 100,000 records are stored in the buffer, with the oldest record being removed when a new record is added.
//! If a PAF record is not found in the buffer, the file is rolled along until the record is found.
use crate::readfish_io::{reader, ByteCounter, DynResult};
use linked_hash_map::LinkedHashMap;
// use rayon::prelude::*;
use std::io::Lines;
use std::{
    io::{BufRead, Read},
    path::{Path, PathBuf},
};
/// Data structure representing sequencing summary information.
///
/// The `SeqSum` struct stores various sequencing summary related fields:
/// - `sequencing_summary_path`: Path to the sequencing summary file.
/// - `writers`: A vector of multiple writers, one for each demultiplexed file.
/// - `record_buffer`: A linked hash map storing the sequencing summary records, with read ID as the key and tuples containing `SeqSumInfo` variants as the values.
/// - `has_barcode`: A boolean indicating whether barcode arrangement is present in the sequencing summary file.
/// - `current_position`: The current position in the file read by the `BufReader`.
/// - `column_indices`: A tuple representing the column indices of `read_id`, `channel`, and `barcode_arrangement` in the sequencing summary file.
///
/// # Examples
/// ```rust,ignore
/// use crate::sequencing_summary::{SeqSum, SeqSumInfo};
/// use std::path::PathBuf;
/// use std::io::Write;
/// use linked_hash_map::LinkedHashMap;
///
/// // Create a new `SeqSum` instance
/// let sequencing_summary_path = PathBuf::from("sequencing_summary.txt");
/// let writers: Vec<Box<dyn Write>> = Vec::new();
/// let record_buffer: LinkedHashMap<String, (SeqSumInfo, SeqSumInfo, SeqSumInfo)> = LinkedHashMap::new();
/// let has_barcode = false;
/// let current_position = 0;
/// let column_indices = (0, 1, 2);
/// let seq_sum = SeqSum {
///     sequencing_summary_path,
///     writers,
///     record_buffer,
///     has_barcode,
///     current_position,
///     column_indices,
/// };
/// ```
pub struct SeqSum {
    /// Path to the sequencing summary file.
    pub sequencing_summary_path: PathBuf,
    /// Multiple writes, one for each demultiplexed file.
    // pub writers: Vec<Box<dyn Write>>,
    /// Record buffer for the sequencing summary
    pub record_buffer: LinkedHashMap<String, (SeqSumInfo, SeqSumInfo, SeqSumInfo)>,
    /// Is barcode_arrangement in this sequencing summary file?
    pub has_barcode: bool,
    /// Current position in file from BufReader
    pub current_position: usize,
    /// Column_indices: (read_id, channel, barcode_arrangement)
    pub column_indices: (usize, usize, usize),
    /// Previous read id. Used to check that we have consumed all of a multiple mapping.
    pub previous_read_id: String,
}

/// Enumeration representing sequenced summary information.
///
/// The `SeqSumInfo` enum holds three possible variants, each representing a different filled:
/// 1. `Channel(usize)`: Stores the channel number of the sequence.
/// 2. `Barcode(String)`: Stores the barcode associated with the sequence.
/// 3. `ReadId(String)`: Stores the unique identifier of the sequence (read ID).
///
/// # Examples
/// ```rust,ignore
/// use crate::SeqSumInfo;
///
/// let channel_info = SeqSumInfo::Channel(3);
/// let barcode_info = SeqSumInfo::Barcode("barcode01".to_string());
/// let read_id_info = SeqSumInfo::ReadId("read12345".to_string());
/// ```
#[derive(Debug, Clone)]
pub enum SeqSumInfo {
    /// Represents a channel with the given usize value.
    Channel(usize),
    /// Represents a barcode with the given String value.
    Barcode(String),
    /// Represents a read ID with the given String value.
    ReadId(String),
}

impl SeqSumInfo {
    /// Get the channel value if the enum variant is Channel, otherwise return None.
    pub fn get_channel(&self) -> Option<usize> {
        if let SeqSumInfo::Channel(channel) = self {
            Some(*channel)
        } else {
            None
        }
    }

    /// Get the barcode value if the enum variant is Barcode, otherwise return None.
    pub fn get_barcode(&self) -> Option<&String> {
        if let SeqSumInfo::Barcode(barcode) = self {
            Some(barcode)
        } else {
            None
        }
    }

    /// Get the read ID value if the enum variant is ReadId, otherwise return None.
    pub fn get_read_id(&self) -> Option<&String> {
        if let SeqSumInfo::ReadId(read_id) = self {
            Some(read_id)
        } else {
            None
        }
    }
}

impl SeqSum {
    /// Create a `SeqSum` instance from a sequencing summary file.
    ///
    /// This function takes the path to a sequencing summary file (`sequencing_summary_path`)
    /// and constructs a [`SeqSum`] instance by parsing the file.
    /// The `SeqSum` struct stores information about the sequencing summary file, including the path, multiple writers for demultiplexed files,
    /// a record buffer, and other relevant data.
    ///
    /// # Arguments
    ///
    /// * `sequencing_summary_path`: An implementation of the `AsRef<Path>` trait that represents the path to the sequencing summary file.
    ///  It can be either a string or a `PathBuf`.
    ///
    /// # Errors
    ///
    /// This function returns a `DynResult<SeqSum>`, which is a type alias for `Result<SeqSum, Box<dyn Error + 'static>>`.
    /// It can return an error if there are issues reading or parsing the sequencing summary file.
    ///
    /// # Examples
    /// ```rust,ignore
    /// # use my_crate::{SeqSum, SeqSumInfo};
    /// # use my_crate::DynResult;
    /// #
    /// # fn main() -> DynResult<()> {
    /// let sequencing_summary_path = "path/to/sequencing_summary.tsv";
    /// let seq_sum = SeqSum::from_file(sequencing_summary_path)?;
    ///
    /// // Now you can use `seq_sum` to perform various operations on the sequencing summary data.
    /// #
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_file(sequencing_summary_path: impl AsRef<Path>) -> DynResult<SeqSum> {
        let sequencing_summary_path = sequencing_summary_path.as_ref().to_path_buf();
        // let writers = vec![];

        let reader = reader(&sequencing_summary_path, None);
        let mut reader = ByteCounter::new(reader);
        let mut lines: Lines<&mut ByteCounter<Box<dyn BufRead + Send>>> = reader.by_ref().lines();
        let headers = lines.next();
        // TODO rewrite function to get the index of a column header
        let read_id_index = headers
            .as_ref()
            .unwrap()
            .as_ref()
            .unwrap()
            .split('\t')
            .position(|column_header| column_header == "read_id");
        let barcode_index = headers
            .as_ref()
            .unwrap()
            .as_ref()
            .unwrap()
            .split('\t')
            .position(|column_header| column_header == "barcode_arrangement");
        let channel_index = headers
            .as_ref()
            .unwrap()
            .as_ref()
            .unwrap()
            .split('\t')
            .position(|column_header| column_header == "channel");
        assert!(
            read_id_index.is_some() && channel_index.is_some(),
            "read_id column header not found in sequencing summary. Header row is likely missing from sequencing summary file."
        );
        let lines_iter = lines.take(100000);
        let processed_lines = LinkedHashMap::from_iter(lines_iter.map(|line| {
            if let Ok(line_content) = line {
                // Process the line content here
                let key = line_content
                    .split('\t')
                    .nth(read_id_index.unwrap())
                    .unwrap()
                    .to_string();
                let selected_elements: Vec<_> = line_content
                    .split('\t')
                    .enumerate()
                    .filter(|(index, _)| {
                        [
                            read_id_index.unwrap_or(usize::MAX),
                            channel_index.unwrap_or(usize::MAX),
                            barcode_index.unwrap_or(usize::MAX),
                        ]
                        .contains(index)
                    })
                    .map(|(_, value)| value)
                    .collect();
                (
                    key,
                    (
                        SeqSumInfo::ReadId(selected_elements[0].to_string()),
                        SeqSumInfo::Channel(selected_elements[1].parse().unwrap()),
                        SeqSumInfo::Barcode(
                            selected_elements
                                .get(2)
                                .unwrap_or(&"no_barcode")
                                .to_string(),
                        ),
                    ),
                )
            } else {
                // Handle any errors that occurred while reading the line
                panic!("failed to read sequencing summary line");
            }
        }));

        Ok(SeqSum {
            sequencing_summary_path,
            // writers,
            record_buffer: processed_lines,
            has_barcode: barcode_index.is_some(),
            current_position: reader.bytes_read(),
            column_indices: (
                read_id_index.unwrap(),
                channel_index.unwrap(),
                barcode_index.unwrap_or(usize::MAX),
            ),
            previous_read_id: String::new(),
        })
    }
    /// Roll along the sequencing summary file until a specific record with the given Read ID is found.
    ///
    /// This function reads the sequencing summary file starting from the current position and searches for
    /// a specific record with the provided Read ID (`query_record_read_id`). If the record is found,
    /// the function stops rolling and the current position in the file is updated.
    ///
    /// The `record_buffer` of the `SeqSum` struct is used to store the sequencing summary records as a linked hash map,
    /// with the Read ID as the key and tuples containing `SeqSumInfo` variants as the values.
    /// The buffer holds a maximum of 100,000 records, and the oldest record is removed when a new record is added.
    ///
    /// # Arguments
    ///
    /// * `query_record_read_id`: A `String` representing the Read ID of the record to search for in the sequencing summary file.
    ///
    /// # Errors
    ///
    /// This function returns a `DynResult<()>` which is a type alias for `Result<(), Box<dyn Error + 'static>>`.
    /// It can return an error if there is an issue reading the sequencing summary file.
    ///
    /// # Examples
    /// ```rust,ignore
    /// # use my_crate::SeqSum;
    /// # use my_crate::DynResult;
    /// #
    /// # fn main() -> DynResult<()> {
    /// let mut seq_sum = SeqSum::from_file("sequencing_summary.tsv")?;
    /// seq_sum.roll_along_file("read123")?;
    /// # Ok(())
    /// # }
    /// ```
    fn roll_along_file(&mut self, query_record_read_id: String) -> DynResult<()> {
        let mut reader = ByteCounter::new(reader(
            &self.sequencing_summary_path,
            Some(self.current_position),
        ));
        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            // do something with line
            let selected_elements: Vec<_> = line
                .split('\t')
                .enumerate()
                .filter(|(index, _)| {
                    [
                        self.column_indices.0,
                        self.column_indices.1,
                        self.column_indices.2,
                    ]
                    .contains(index)
                })
                .map(|(_, value)| value)
                .collect();

            let key = line
                .split('\t')
                .nth(self.column_indices.0)
                .unwrap()
                .to_string();
            self.record_buffer.pop_front().unwrap();
            self.record_buffer.insert(
                key,
                (
                    SeqSumInfo::ReadId(selected_elements[0].to_string()),
                    SeqSumInfo::Channel(selected_elements[1].parse().unwrap()),
                    SeqSumInfo::Barcode(
                        selected_elements
                            .get(2)
                            .unwrap_or(&"no_barcode")
                            .to_string(),
                    ),
                ),
            );
            if *selected_elements[0] == query_record_read_id {
                break;
            }
            line.clear();
        }
        self.current_position += reader.bytes_read();
        Ok(())
    }

    /// Get the sequencing summary record associated with the given `query_name`.
    /// The record is returned as a tuple containing three `SeqSumInfo` variants.
    ///
    /// This function searches for the sequencing summary record corresponding to the provided `query_name` in the record buffer of the `SeqSum` struct.
    ///  If the record is found in the buffer, it is returned.
    /// If not, the function rolls along the sequencing summary file to find the record with the matching Read ID (`query_name`).
    ///
    /// # Arguments
    ///
    /// * `query_name`: A `&str` representing the Read ID to search for in the sequencing summary records.
    /// * `previous_read_id`: An optional mutable reference to a `String` representing the Read ID of the previous record.
    /// This is used to keep track of the current position in the sequencing summary file.
    ///
    /// # Errors
    ///
    /// This function returns a `DynResult<(SeqSumInfo, SeqSumInfo, SeqSumInfo)>`,
    /// which is a type alias for `Result<(SeqSumInfo, SeqSumInfo, SeqSumInfo), Box<dyn Error + 'static>>`.
    /// It can return an error if there is an issue reading the sequencing summary file while rolling along
    /// to find the record with the matching Read ID.
    ///
    /// # Examples
    /// ```rust,ignore
    /// # use my_crate::{SeqSum, SeqSumInfo};
    /// # use my_crate::DynResult;
    /// #
    /// # fn main() -> DynResult<()> {
    /// let mut seq_sum = SeqSum::new("sequencing_summary.tsv")?;
    ///
    /// let query_name = "read123";
    /// let previous_read_id = Some(&mut String::new());
    ///
    /// let record = seq_sum.get_record(query_name, previous_read_id)?;
    /// match record {
    ///     (SeqSumInfo::ReadId(read_id), SeqSumInfo::Channel(channel), SeqSumInfo::Barcode(barcode)) => {
    ///         println!("Read ID: {}", read_id);
    ///         println!("Channel: {}", channel);
    ///         println!("Barcode: {}", barcode);
    ///     },
    ///     _ => println!("Record not found."),
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_record(
        &mut self,
        query_name: &str,
        previous_query_name: Option<&str>,
    ) -> DynResult<(SeqSumInfo, SeqSumInfo, SeqSumInfo)> {
        if (query_name != previous_query_name.unwrap_or(&self.previous_read_id))
            & (!previous_query_name.unwrap_or("").is_empty())
        {
            self.record_buffer
                .remove(previous_query_name.unwrap())
                .unwrap();
        } else {
        }
        match self.record_buffer.get(query_name) {
            Some(record) => Ok(record.clone()),
            None => {
                // Assuming multiple mappings are in a block in a PAF file
                self.roll_along_file(query_name.to_string())?;
                Ok(self.record_buffer.get(query_name).unwrap().clone())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_resource_dir() -> PathBuf {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("resources");
        path
    }

    fn get_test_file(file: &str) -> PathBuf {
        let mut path = get_resource_dir();
        path.push(file);
        path
    }

    #[test]
    fn test_seq_sum_from_file() {
        let seq_sum_file_path = get_test_file("seq_sum_PAK09329.txt");
        let seq_sum = SeqSum::from_file(seq_sum_file_path).unwrap();
        assert_eq!(
            seq_sum.sequencing_summary_path,
            get_test_file("seq_sum_PAK09329.txt")
        );
        assert_eq!(seq_sum.record_buffer.len(), 100000);
        assert!(seq_sum.has_barcode);
    }
}
