//! Paf file functions
//! In this module we implement a Paf struct and functions to read and write Paf files.
//! A lot of this was lifted from https://github.com/mrvollger/rustybam/blob/main/src/paf.rs
//!

use crate::{
    readfish::Conf,
    readfish_io::{reader, DynResult},
    sequencing_summary::SeqSum,
    Summary,
};
use lazy_static::lazy_static;
use regex::Regex;
use std::{
    io::BufRead,
    path::{Path, PathBuf},
};

lazy_static! {
    static ref PAF_TAG: Regex = Regex::new("(..):(.):(.*)").unwrap();
}

/// Store metadata that is provided by a tuple in a call to parse_paf_by_iter in lib.rs.
/// See also `[sequencing_summary::SeqSumInfo]`.
#[derive(Debug)]
pub struct Metadata {
    /// The identifier for the read.
    pub read_id: String,
    /// The channel associated with the read.
    pub channel: usize,
    /// An optional barcode associated with the read, if available.
    pub barcode: Option<String>,
}

impl From<(String, usize, Option<String>)> for Metadata {
    fn from(value: (String, usize, Option<String>)) -> Self {
        Metadata {
            read_id: value.0,
            channel: value.1,
            barcode: value.2,
        }
    }
}

impl Metadata {
    /// Get the identifier for the read.
    pub fn read_id(&self) -> &String {
        &self.read_id
    }

    /// Get the channel associated with the read.
    pub fn channel(&self) -> usize {
        self.channel
    }

    /// Get the optional barcode associated with the read, if available.
    pub fn barcode(&self) -> Option<&String> {
        self.barcode.as_ref()
    }
}

/// Store a PafRecord for quick unpacking to update the summary
#[derive(Debug, Clone)]
pub struct PafRecord {
    /// The name of the query sequence (read).
    pub query_name: String,
    /// The length of the query sequence (read).
    pub query_length: usize,
    /// The start position of the alignment on the query sequence (read).
    pub query_start: usize,
    /// The end position of the alignment on the query sequence (read).
    pub query_end: usize,
    /// The strand of the alignment ('+' or '-').
    pub strand: char,
    /// The name of the target sequence (reference).
    pub target_name: String,
    /// The length of the target sequence (reference).
    pub target_length: usize,
    /// The start position of the alignment on the target sequence (reference).
    pub target_start: usize,
    /// The end position of the alignment on the target sequence (reference).
    pub target_end: usize,
    /// The number of matching bases in the alignment.
    pub nmatch: usize,
    /// The total length of the alignment.
    pub aln_len: usize,
    /// The mapping quality of the alignment.
    pub mapq: usize,
    // pub cigar: CigarString,
    // A vector of additional tags associated with the alignment.
    // pub tags: Vec<String>,
    // pub tpos_aln: Vec<u64>,
    // pub qpos_aln: Vec<u64>,
    // pub long_cigar: CigarString,
    // pub id: String,
    // pub order: u64,
    // pub contained: bool,
}
/// Errors that can occur while parsing PAF (Pairwise mApping Format) files.
#[derive(Debug)]
pub enum Error {
    /// An error occurred while parsing the CIGAR string in the PAF record.
    PafParseCigar {
        /// The error message.
        msg: String,
    },
    /// An error occurred while parsing the CS (Coordinate System) tag in the PAF record.
    PafParseCS {
        /// The error message.
        msg: String,
    },
    /// An error occurred while parsing an integer value in the PAF record.
    ParseIntError {
        /// The error message.
        msg: String,
    },
    /// An error occurred while parsing a column in the PAF record.
    ParsePafColumn {},
}

/// A type alias for a Result with the error type specialized to `crate::paf::Error`.
type PafResult<T> = Result<T, crate::paf::Error>;

impl PafRecord {
    /// New paf record
    pub fn new(t: Vec<&str>) -> PafResult<PafRecord> {
        // make the record
        let rec = PafRecord {
            query_name: t[0].to_string(),
            query_length: t[1]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            query_start: t[2]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            query_end: t[3]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            strand: t[4].parse::<char>().map_err(|_| Error::ParsePafColumn {})?,
            target_name: t[5].to_string(),
            target_length: t[6]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            target_start: t[7]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            target_end: t[8]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            nmatch: t[9]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            aln_len: t[10]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
            mapq: t[11]
                .parse::<usize>()
                .map_err(|_| Error::ParsePafColumn {})?,
        };
        Ok(rec)
    }
}

/// A struct representing a PAF record reader and writers for demultiplexing.
///
/// This struct holds a reader and a list of writers used for demultiplexing PAF records
/// into different files. The `reader` field is a `Box<dyn BufRead + Send>` representing a
/// buffered input reader from which PAF records are read. The `writers` field is a `Vec<Box<dyn Write>>`
/// holding multiple output writers for writing the demultiplexed PAF records into different files.
///
/// # Fields
///
/// * `reader`: A boxed trait object implementing `BufRead` and `Send`, used as the input reader
///   for reading PAF records.
/// * `writers`: A vector of boxed trait objects implementing `Write`, used as the output writers
///   for writing the demultiplexed PAF records into different files.
/// * `paf_file`: The path to the PAF file.
///
/// # Examples
///
/// ```rust, ignore
/// use std::fs::File;
/// use std::io::{BufReader, BufWriter};
/// use std::path::Path;
///
/// // Create a reader for the PAF file
/// let file_path = Path::new("example.paf");
/// let file = File::open(file_path).expect("Error: Failed to open file");
/// let reader = Box::new(BufReader::new(file));
///
/// // Create multiple writers for demultiplexing the PAF records
/// let writer1 = Box::new(BufWriter::new(File::create("output1.paf").unwrap()));
/// let writer2 = Box::new(BufWriter::new(File::create("output2.paf").unwrap()));
/// let writers = vec![writer1, writer2];
///
/// // Create a PAF object
/// let paf = Paf { reader, writers };
/// ```
///
pub struct Paf {
    /// The provided PAF file.
    pub paf_file: PathBuf,
    /// Reader for the Paf file.
    pub reader: Box<dyn BufRead + Send>,
    // / Multiple writes, one for each demultiplexed file.
    // pub writers: Vec<Box<dyn Write>>,
}

impl Paf {
    /// Create a new `Paf` object with the given PAF file.
    ///
    /// This function creates a new `Paf` object by parsing the specified PAF file
    /// and initializing the `reader` field with the resulting buffered input reader.
    /// The `writers` field is initialized as an empty vector of output writers.
    ///
    /// # Arguments
    ///
    /// * `paf_file`: An implementation of the `AsRef<Path>` trait representing the path to the PAF file.
    ///
    /// # Returns
    ///
    /// A new `Paf` object with the parsed PAF file as the input reader and an empty vector of writers.
    ///
    /// # Panics
    ///
    /// This function will panic if there is an error while parsing the PAF file or creating the buffered input reader.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use std::path::Path;
    /// use readfish_tools::Paf;
    ///
    /// // Create a new Paf object from the "example.paf" file
    /// let paf_file_path = Path::new("example.paf");
    /// let paf = Paf::new(paf_file_path);
    /// ```
    ///
    pub fn new(paf_file: impl AsRef<Path>) -> Paf {
        Paf {
            paf_file: paf_file.as_ref().to_path_buf(),
            reader: open_paf_for_reading(paf_file).unwrap(),
            // writers: vec![],
        }
    }

    /// Demultiplexes the PAF file by processing each line and obtaining corresponding sequencing summary records.
    ///
    /// This function reads the PAF file line by line, parses each line, and processes the custom tags present in the PAF format.
    /// These custom tags are add by readfish's implementation summarise on the Aligner.
    /// If the `sequencing_summary` argument is provided, it retrieves the sequencing summary record for each line's query name.
    /// The function processes custom tags in the PAF file and ensures they are present. If `sequencing_summary` is None and custom tags are missing,
    /// the function will panic.
    ///
    /// If `sequencing_summary` is provided, the function retrieves the sequencing summary record for each query name using the `get_record` function.
    /// If a sequencing summary record is not found in the buffer, the function reads from the sequencing summary file until the record is found.
    /// The function consumes the bytes in the PAF file and updates the `previous_read_id` to avoid removing multiple mappings from the `sequencing_summary`
    /// only when the new Read Id is not the same as the old read_id.
    ///
    /// # Arguments
    ///
    /// - `toml`: A reference to the `Conf` struct, which contains configuration settings.
    /// - `sequencing_summary`: An optional mutable reference to the `SeqSum` struct, representing the sequencing summary file.
    ///
    /// # Errors
    ///
    /// This function returns a `DynResult`, which is a specialized `Result` type with an error message.
    /// An error is returned if there is any issue reading the PAF file or if the sequencing summary file is not found.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// // Import necessary libraries
    /// use std::error::Error;
    /// use my_crate::{SeqSum, Conf};
    ///
    /// // Create a new sequencing summary instance
    /// let mut sequencing_summary = SeqSum::from_file("path/to/sequencing_summary.toml")?;
    ///
    /// // Load the TOML configuration
    /// let toml = Conf::from_file("path/to/config.toml")?;
    ///
    /// // Demultiplex the PAF file using the sequencing summary
    /// sequencing_summary.demultiplex(&toml, Some(&mut sequencing_summary))?;
    /// ```
    pub fn demultiplex(
        &mut self,
        _toml: &mut Conf,
        sequencing_summary: Option<&mut SeqSum>,
        mut summary: Option<&mut Summary>,
    ) -> DynResult<()> {
        let seq_sum = sequencing_summary.unwrap();

        // Remove multiple mappings from seq_sum dictionary only when the new Read Id is not the same as the old read_id
        for line in open_paf_for_reading(self.paf_file.clone())?.lines() {
            let (paf_record, read_on, condition_name) =
                _parse_paf_line(line?, _toml, None, Some(seq_sum))?;

            if let Some(summary) = summary.as_deref_mut() {
                let condition_summary = summary.conditions(condition_name.as_str());
                condition_summary.update(paf_record, read_on).unwrap();
            }
        }
        Ok(())
    }
}

/// Parses the PAF file and returns a buffered reader for further processing.
///
/// This function takes the `file_name` as an input and returns a `Result` containing
/// a boxed buffered reader (`Box<dyn BufRead + Send>`) if the PAF file is valid.
///
/// # Arguments
///
/// * `file_name`: The path to the PAF file to be parsed. It should implement the `AsRef<Path>` trait.
///
/// # Returns
///
/// A `Result` containing the buffered reader (`Box<dyn BufRead + Send>`) if the PAF file is valid,
/// otherwise an `Err` containing an error message.
///
/// # Errors
///
/// The function returns an `Err` in the following cases:
/// - If the file is empty (contains no bytes), it returns an error with the message "Error: empty file".
/// - If the file format is invalid, specifically if any of the first twelve columns contains a ':' character,
///   it returns an error with the message "Error: invalid format for PAF file. Missing one of first twelve columns, or values contain a :."
/// - If there are any I/O errors while reading the file, the function returns an error with the specific I/O error message.
///
/// # Example
///
/// ```rust,ignore
/// use std::path::Path;
/// use std::io::BufRead;
///
/// let file_name = Path::new("example.paf");
/// match open_paf_for_reading(file_name) {
///     Ok(reader) => {
///         let mut lines = reader.lines();
///         // Process the PAF file line by line
///         for line in lines {
///             match line {
///                 Ok(line_content) => {
///                     // Process the content of the PAF line
///                     // ...
///                 }
///                 Err(error) => {
///                     eprintln!("Error while reading PAF line: {}", error);
///                 }
///             }
///         }
///     }
///     Err(error) => {
///         eprintln!("Error while parsing PAF file: {}", error);
///     }
/// }
/// ```
pub fn open_paf_for_reading(file_name: impl AsRef<Path>) -> DynResult<Box<dyn BufRead + Send>> {
    // create reader to check file first line
    let mut paf_file = reader(&file_name, None);

    // Check the file isn't empty
    let mut buffer = [0; 1];
    let bytes_read = paf_file.read(&mut buffer)?;
    if bytes_read == 0 {
        return Err("Error: empty file".into());
    }
    let mut line = String::new();
    paf_file.read_line(&mut line)?;
    let t: Vec<&str> = line.split_ascii_whitespace().collect();
    if t.iter().take(12).any(|item| item.contains(':')) {
        return Err("Error: invalid format for Paf file. Missing one of first twelve columns, or values contain a :.".into());
    }

    let paf_file = reader(file_name, None);

    Ok(paf_file)
}

/// Parses a line from the PAF file and extracts relevant information.
///
/// This function takes a PAF line (as a reference to a string) and attempts to parse it to extract
/// relevant information, including creating a [`PafRecord`] and making decisions based on the provided
/// metadata or sequencing summary. It returns a tuple containing the `PafRecord`, a boolean value
/// indicating if the read is considered "on-target", and the condition name associated with the read.
///
/// # Arguments
///
/// * `paf_line`: A reference to a string slice representing a line from the PAF file.
/// * `_toml`: A reference to a `Conf` struct, holding configuration information.
/// * `meta_data`: An optional mutable reference to a `Metadata` struct containing read metadata.
/// * `sequencing_summary`: An optional mutable reference to a `SeqSum` struct containing sequencing summary data.
///
/// # Returns
///
/// A `DynResult` holding a tuple containing the following elements:
/// * `PafRecord`: The parsed PAF record representing the alignment information.
/// * `bool`: A boolean value indicating if the read is considered "on-target".
/// * `&'a String`: A reference to the condition name associated with the read.
///
/// # Panics
///
/// This function panics if the PAF line contains missing items in the first 12 columns or if both `meta_data`
/// and `sequencing_summary` are `None`.
///
/// # Examples
///
/// ```rust, ignore
/// # use your_crate::{PafRecord, Metadata, SeqSum, _parse_paf_line};
///
/// // Assuming we have valid inputs
/// let paf_line = "read123 200 0 200 + contig123 300 0 300 200 200 50 ch=1";
/// let _toml = Conf::default();
/// let mut metadata = Metadata {
///     read_id: "read123".to_string(),
///     channel: 1,
///     barcode: Some("sampleA".to_string()),
/// };
/// let mut seq_sum = SeqSum::default();
///
/// let result = _parse_paf_line(paf_line, &_toml, Some(&mut metadata), Some(&mut seq_sum));
///
/// match result {
///     Ok((paf_record, read_on, condition_name)) => {
///         // Do something with the parsed data
///     }
///     Err(err) => {
///         // Handle the error
///     }
/// }
/// ```
pub fn _parse_paf_line<'a>(
    paf_line: impl AsRef<str>,
    _toml: &'a Conf,
    meta_data: Option<&mut Metadata>,
    sequencing_summary: Option<&mut SeqSum>,
) -> DynResult<(PafRecord, bool, &'a String)> {
    let line = paf_line.as_ref();
    let t: Vec<&str> = line.split_ascii_whitespace().collect();
    // Todo do without clone
    let paf_record = PafRecord::new(t.clone()).unwrap();
    // Check first 12 columns for missing items, assumes tags will have been brought forwards
    assert!(
        t.iter().take(12).all(|item| !item.contains(':')),
        "Missing colon in PAF line: {}",
        line
    );
    // check if we have custom tags from readfish aligner analyse
    let channel: usize;
    let barcode: Option<String>;
    // for token in t.iter().skip(12) {
    //     debug_assert!(PAF_TAG.is_match(token));
    //     let caps = PAF_TAG.captures(token).unwrap();
    //     let tag = &caps[1];
    //     // let value = &caps[3];
    //     if (tag == "ch") | (tag == "ba") {
    //         has_tags = true;
    //     }
    // }
    // Break the Paf line into its components
    let query_name = t[0];
    // let query_length: usize = t[1].parse()?;
    let strand = t[4];
    let contig = t[5];
    // let contig_length: usize = t[6].parse()?;
    let mapping_start: usize = t[7].parse()?;
    let read_on: bool;
    if meta_data.is_none() & sequencing_summary.is_none() {
        panic!("Cannot parse paf line without provided metdata or sequencing summary_file");
    }
    // If sequencing summary is provided, get the sequencing summary record for the query name
    // Use it for things like barcodes and channels
    if let Some(seq_sum_struct) = sequencing_summary {
        let seq_sum_record = seq_sum_struct.get_record(query_name, None);
        if let Ok(record) = seq_sum_record {
            read_on = _toml.make_decision(
                record.1.get_channel().unwrap(),
                record.2.get_barcode().map(|x| x.as_str()),
                contig,
                strand,
                mapping_start,
            );
            channel = record.1.get_channel().unwrap();
            barcode = Some(record.2.get_barcode().unwrap_or(&"".to_string()).clone());
        } else {
            return Err("Error: sequencing summary record not found".into());
        }
        seq_sum_struct.previous_read_id = query_name.to_string();
    // We must have metatdata
    } else {
        let metadata = meta_data.unwrap();
        // println!("{contig}, {strand}, {mapping_start}");
        read_on = _toml.make_decision(
            metadata.channel(),
            metadata.barcode().map(|x| x.as_str()),
            contig,
            strand,
            mapping_start,
        );
        channel = metadata.channel();
        barcode = Some(metadata.barcode().unwrap_or(&"".to_string()).clone());
    }
    // get the condition so we can access name etc.
    let (_control, condition) = _toml.get_conditions(channel, barcode)?;
    let condition = condition.get_condition();
    let condition_name = &condition.name;

    Ok((paf_record, read_on, condition_name))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_resource_dir() -> PathBuf {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("resources/");
        path
    }

    fn get_test_file(file: &str) -> PathBuf {
        let mut path = get_resource_dir();
        path.push(file);
        path
    }

    #[test]
    fn test_from_tuple() {
        let tuple = ("ABC123".to_string(), 1, Some("BCDE".to_string()));
        let metadata = Metadata::from(tuple);

        assert_eq!(metadata.read_id(), "ABC123");
        assert_eq!(metadata.channel(), 1);
        assert_eq!(metadata.barcode(), Some(&"BCDE".to_string()));
    }

    #[test]
    fn test_from_tuple_no_barcode() {
        let tuple = ("XYZ789".to_string(), 2, None);
        let metadata = Metadata::from(tuple);

        assert_eq!(metadata.read_id(), "XYZ789");
        assert_eq!(metadata.channel(), 2);
        assert_eq!(metadata.barcode(), None);
    }

    #[test]
    fn test_read_id() {
        let metadata = Metadata {
            read_id: "ABC123".to_string(),
            channel: 1,
            barcode: None,
        };

        assert_eq!(metadata.read_id(), "ABC123");
    }

    #[test]
    fn test_channel() {
        let metadata = Metadata {
            read_id: "ABC123".to_string(),
            channel: 1,
            barcode: Some("BCDE".to_string()),
        };

        assert_eq!(metadata.channel(), 1);
    }

    #[test]
    fn test_barcode_present() {
        let metadata = Metadata {
            read_id: "ABC123".to_string(),
            channel: 1,
            barcode: Some("BCDE".to_string()),
        };

        assert_eq!(metadata.barcode(), Some(&"BCDE".to_string()));
    }

    #[test]
    fn test_barcode_absent() {
        let metadata = Metadata {
            read_id: "ABC123".to_string(),
            channel: 1,
            barcode: None,
        };

        assert_eq!(metadata.barcode(), None);
    }

    #[test]
    fn test_from_file_valid_paf() {
        let file_name = get_test_file("test_hum_4000.paf");
        let result = open_paf_for_reading(file_name);
        assert!(
            result.is_ok(),
            "Expected Ok, but got an error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_from_file_invalid_paf() {
        let file_name = get_test_file("invalid_file.paf");
        let result = open_paf_for_reading(file_name);
        assert!(result.is_err(), "Expected Err, but got Ok");
    }

    #[test]
    fn test_from_file_empty_file() {
        let file_name = get_test_file("empty.paf");
        let result = open_paf_for_reading(file_name);
        assert!(result.is_err(), "Expected Err, but got Ok");
    }

    #[test]
    #[should_panic]
    fn test_from_file_nonexistent_file() {
        let file_name = get_test_file("no_existo.paf");
        let result = open_paf_for_reading(file_name);
        assert!(result.is_err(), "Expected Err, but got Ok");
    }

    #[test]
    fn test_paf_from_file() {
        open_paf_for_reading(get_test_file("test_hum_4000.paf")).unwrap();
        // assert_eq!(paf.records.len(), 4148usize);
    }
}
