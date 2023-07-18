//! Paf file functions
//! In this module we implement a Paf struct and functions to read and write Paf files.
//! A lot of this was lifted from https://github.com/mrvollger/rustybam/blob/main/src/paf.rs
//!
use flate2::{read, Compression};
use gzp::{deflate::Bgzf, ZBuilder};
use lazy_static::lazy_static;
use regex::Regex;
use std::{
    error::Error,
    ffi::OsStr,
    fs::File,
    io::{stdin, stdout, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

lazy_static! {
    static ref PAF_TAG: Regex = Regex::new("(..):(.):(.*)").unwrap();
}

/// Small default BUFFER_SIZE for buffered readers
const BUFFER_SIZE: usize = 32 * 1024;

/// Dynamic result type for holding either a generic value or an error
type DynResult<T> = Result<T, Box<dyn Error + 'static>>;

/// Get a buffered input reader from stdin or a file
///
/// # Arguments
///
/// * `path`: Optional path to a file. If provided, the function will attempt to open the file and create a buffered reader from it. If set to `None`, the function will create a buffered reader from stdin.
///
/// # Returns
///
/// A dynamic result that represents either a buffered reader on success or an error wrapped in a `Box<dyn Error + 'static>` on failure.
///
/// # Examples
///
/// ```rust,ignore
/// let reader = _get_reader_from_path(Some(PathBuf::from("path/to/file")))?;
/// ```
fn _get_reader_from_path(path: Option<PathBuf>) -> DynResult<Box<dyn BufRead + Send + 'static>> {
    let reader: Box<dyn BufRead + Send + 'static> = match path {
        Some(path) => {
            // stdin
            if path.as_os_str() == "-" {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin()))
            } else {
                // open file
                Box::new(BufReader::with_capacity(BUFFER_SIZE, File::open(path)?))
            }
        }
        // nothing passed as input, read from stdin
        None => Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin())),
    };
    Ok(reader)
}

/// Read normal or compressed files seamlessly
///
/// This function provides a convenient way to read both normal and compressed files.
///  It automatically detects whether the file is compressed based on the presence of a
/// `.gz` or `.bgz` extension in the filename.
///
/// # Examples
///
/// Reading from an uncompressed file:
///
/// ```rust,ignore
/// use std::io::BufRead;
///
/// let n_lines = reader("file.txt").lines().count();
/// assert_eq!(n_lines, 10);
/// ```
///
/// Reading from a compressed file:
///
/// ```rust,ignore
/// use std::io::BufRead;
///
/// let n_lines_gz = reader("file.txt.gz").lines().count();
/// assert_eq!(n_lines_gz, 10);
/// ```
///
/// In the examples above, the `reader` function seamlessly handles both uncompressed and compressed files.
///  It returns a buffered reader (`Box<dyn BufRead + Send + 'static>`) that can be used to read the file's contents line by line.
///
/// # Arguments
///
/// * `filename`: The path or filename of the file to read. If "-" is provided, the function will read from stdin.
///
/// # Returns
///
/// A boxed trait object implementing `BufRead`, which can be used to read the contents of the file.
/// Uses the presence of a `.gz` or `.bgz` extension to decide
pub fn reader(filename: impl AsRef<Path>) -> Box<dyn BufRead + Send + 'static> {
    let ext = filename.as_ref().extension();
    let path: PathBuf = filename.as_ref().to_path_buf();
    // Handle Gzipped files first, since need to use flate2::read::GzDecoder
    if ext == Some(OsStr::new("gz")) {
        let file = match File::open(&path) {
            Err(why) => panic!("couldn't open {}: {}", path.display(), why),
            Ok(file) => file,
        };
        Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            read::GzDecoder::new(file),
        ))
    // } else if ext == Some(OsStr::new("bgz")) {
    //     Box::new(BufReader::new(BgzfSyncReader::new(
    //         get_input(Some(path)).expect("Error: cannot read input file."),
    //     )))
    } else {
        _get_reader_from_path(Some(path)).expect("Error: cannot read input file")
    }
}

/// Gets a buffered output writer from stdout or a file.
///
/// This function creates a buffered output writer from either stdout or a file specified
/// by the provided `path`. If the `path` is [`Some`], it creates a buffered writer for the
/// specified file. If the `path` is `None`, it creates a buffered writer for stdout.
///
/// The function returns a [`Result`] containing the boxed writer if successful, or an error
/// if the file cannot be created or if an I/O error occurs.
///
/// # Arguments
///
/// * `path` - An optional path to the file. If `Some`, a buffered writer for the file will be created.
///            If `None`, a buffered writer for stdout will be created.
///
/// # Returns
///
/// A `Result` containing the boxed writer if successful, or an error message if the file cannot be created
/// or if an I/O error occurs.
///
/// # Examples
///
/// ```rust,ignore
/// use std::path::PathBuf;
///
/// let path = Some(PathBuf::from("output.txt"));
/// let writer = get_output(path);
///
/// match writer {
///     Ok(w) => {
///         // Write data using the writer
///     }
///     Err(err) => {
///         eprintln!("Error creating output writer: {}", err);
///     }
/// }
/// ```
fn _get_writer_from_path(
    path: Option<PathBuf>,
) -> Result<Box<dyn Write + Send + 'static>, std::io::Error> {
    let writer: Box<dyn Write + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, stdout()))
            } else {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, File::create(path)?))
            }
        }
        None => Box::new(BufWriter::with_capacity(BUFFER_SIZE, stdout())),
    };
    Ok(writer)
}

/// Write data to normal or compressed files seamlessly.
/// The function determines the file type based on the presence of the `.gz` extension.
///
/// # Arguments
///
/// * `filename` - The name of the file to write to, including the extension.
///
/// # Returns
///
/// A boxed trait object (`Box<dyn Write>`) representing the writer for the specified file.
///
/// # Examples
///
/// ```rust,ignore
/// use std::io::Write;
/// let mut writer = writer("output.txt");
/// writer.write_all(b"Hello, world!").expect("Failed to write data");
/// ```
pub fn writer(filename: &str) -> Box<dyn Write> {
    let ext = Path::new(filename).extension();
    let path = PathBuf::from(filename);
    let buffer = _get_writer_from_path(Some(path)).expect("Error: cannot create output file");

    if ext == Some(OsStr::new("gz")) {
        let writer = ZBuilder::<Bgzf, _>::new()
            .num_threads(8)
            .compression_level(Compression::new(6))
            .from_writer(buffer);
        Box::new(writer)
    } else {
        buffer
    }
}

/// Paf record
pub fn from_file(file_name: impl AsRef<Path>) {
    let paf_file = reader(file_name);
    // let mut paf = Paf::new();
    // read the paf recs into a vector
    for (_index, line) in paf_file.lines().enumerate() {
        println!("Line: {}", line.unwrap());
        // log::trace!("{:?}", line);
        // match PafRecord::new(&line.unwrap()) {
        //     Ok(mut rec) => {
        //         rec.check_integrity().unwrap();
        //         paf.records.push(rec);
        //     }
        //     Err(_) => eprintln!("\nUnable to parse PAF record. Skipping line {}", index + 1),
        // }
        // log::debug!("Read PAF record number: {}", index + 1);
    }
    // paf
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
    fn test_paf_from_file() {
        from_file(get_test_file("test_hum_4000.paf"));
        // assert_eq!(paf.records.len(), 4148usize);
    }

    #[test]
    fn test_reader() {
        let n_lines = reader(get_test_file("test_hum_4000.paf")).lines().count();
        assert_eq!(n_lines, 4148usize);

        let n_lines_gz = reader(get_test_file("test_hum_4000.paf.gz"))
            .lines()
            .count();
        assert_eq!(n_lines_gz, 4148usize);
    }
}
