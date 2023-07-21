//! Paf file functions
//! In this module we implement a Paf struct and functions to read and write Paf files.
//! A lot of this was lifted from https://github.com/mrvollger/rustybam/blob/main/src/paf.rs
//!

use crate::{
    readfish::Conf,
    readfish_io::{reader, DynResult},
    sequencing_summary::SeqSum,
};
use lazy_static::lazy_static;
use regex::Regex;
use std::{
    io::{BufRead, Write},
    path::{Path, PathBuf},
};

lazy_static! {
    static ref PAF_TAG: Regex = Regex::new("(..):(.):(.*)").unwrap();
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
    /// Multiple writes, one for each demultiplexed file.
    pub writers: Vec<Box<dyn Write>>,
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
    /// ```rust
    /// use std::path::Path;
    /// use readfish::Paf;
    ///
    /// // Create a new Paf object from the "example.paf" file
    /// let paf_file_path = Path::new("example.paf");
    /// let paf = Paf::new(paf_file_path);
    /// ```
    ///
    pub fn new(paf_file: impl AsRef<Path>) -> Paf {
        Paf {
            paf_file: paf_file.as_ref().to_path_buf(),
            reader: parse_paf_file(paf_file).unwrap(),
            writers: vec![],
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
        _toml: &Conf,
        mut sequencing_summary: Option<&mut SeqSum>,
    ) -> DynResult<()> {
        // Remove multiple mappings from seq_sum dictionary only when the new Read Id is not the same as the old read_id
        let mut previous_read_id = String::new();
        for (_index, line) in parse_paf_file(self.paf_file.clone())?.lines().enumerate() {
            let line = line?;
            println!("line: {}", line);
            let t: Vec<&str> = line.split_ascii_whitespace().collect();
            assert!(
                t.iter().take(12).all(|item| !item.contains(':')),
                "Missing colon in PAF line: {}",
                line
            );
            println!("t: {:?}", t);
            let mut has_tags: bool = sequencing_summary.is_some();
            for token in t.iter().skip(12) {
                debug_assert!(PAF_TAG.is_match(token));
                let caps = PAF_TAG.captures(token).unwrap();
                let tag = &caps[1];
                // let value = &caps[3];
                if (tag == "ch") | (tag == "ba") {
                    has_tags = true;
                }
            }
            let query_name = t[0];

            // Panic if we don't have our custom tags and the sequencing summary file is None
            if !has_tags & sequencing_summary.is_none() {
                panic!("Missing custom tags in PAF line: {}", line);
            }
            if sequencing_summary.is_some() {
                let seq_sum_struct = sequencing_summary.as_deref_mut().unwrap();
                let seq_sum_record =
                    seq_sum_struct.get_record(query_name, Some(&mut previous_read_id));
                println!(
                    "seq_sum_record: {:?}, query_name: {:#?}",
                    seq_sum_record, query_name
                );
            }
            if previous_read_id.is_empty() {
                previous_read_id = query_name.to_string();
            }
        }
        Ok(())
    }
}

/// Reads and parses a PAF file, extracting relevant information from each line.
///
/// A PAF (Pairwise mApping Format) file contains tab-separated columns, each representing
/// alignment information between sequences. The function reads the file line by line,
/// parses each line, and collects the necessary information for further processing.
///
/// # Arguments
///
/// * `file_name`: An implementation of the `AsRef<Path>` trait, representing the path to the PAF file.
///
/// # Errors
///
/// The function returns a result containing either `Ok(())` if the PAF file is successfully read
/// and parsed, or an error message in case of any issues with file access
/// or incorrect PAF format.
///
/// # Examples
///
/// ```rust,no_run
/// use std::path::Path;
/// use my_crate::from_file;
///
/// // Provide the path to the PAF file
/// let file_path = Path::new("path/to/your_paf_file.paf");
///
/// // Call the function to parse the PAF file
/// let result = from_file(file_path);
///
/// // Check the result
/// match result {
///     Ok(_) => println!("PAF file parsed successfully!"),
///     Err(err) => println!("Error: {}", err),
/// }
/// ```
pub fn parse_paf_file(file_name: impl AsRef<Path>) -> DynResult<Box<dyn BufRead + Send>> {
    let mut paf_file = reader(&file_name, None);

    // Check the file isn't empty
    let mut buffer = [0; 1];
    let bytes_read = paf_file.read(&mut buffer)?;
    let paf_file = reader(file_name, None);
    if bytes_read == 0 {
        return Err("Error: empty file".into());
    }
    Ok(paf_file)
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
    fn test_from_file_valid_paf() {
        let file_name = get_test_file("test_hum_4000.paf");
        let result = parse_paf_file(file_name);
        assert!(
            result.is_ok(),
            "Expected Ok, but got an error: {:?}",
            result.err()
        );
    }

    #[test]
    fn test_from_file_invalid_paf() {
        let file_name = get_test_file("invalid_file.paf");
        let result = parse_paf_file(file_name);
        assert!(result.is_err(), "Expected Err, but got Ok");
    }

    #[test]
    fn test_from_file_empty_file() {
        let file_name = get_test_file("empty.paf");
        let result = parse_paf_file(file_name);
        assert!(result.is_err(), "Expected Err, but got Ok");
    }

    #[test]
    #[should_panic]
    fn test_from_file_nonexistent_file() {
        let file_name = get_test_file("no_existo.paf");
        let result = parse_paf_file(file_name);
        assert!(result.is_err(), "Expected Err, but got Ok");
    }

    #[test]
    fn test_paf_from_file() {
        parse_paf_file(get_test_file("test_hum_4000.paf")).unwrap();
        // assert_eq!(paf.records.len(), 4148usize);
    }
}
