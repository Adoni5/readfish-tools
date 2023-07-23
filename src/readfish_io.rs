//! Io functions for reading/writing gzipped or uncompressed files.

use flate2::{read::GzDecoder, Compression};
use gzp::{deflate::Bgzf, ZBuilder};
use std::{
    error::Error,
    ffi::OsStr,
    fs::File,
    io::{
        stdin, stdout, BufRead, BufReader, BufWriter, Read, Result as ioResult, Seek, SeekFrom,
        Write,
    },
    path::{Path, PathBuf},
};
/// Small default BUFFER_SIZE for buffered readers
const BUFFER_SIZE: usize = 32 * 1024;

/// Dynamic result type for holding either a generic value or an error
pub type DynResult<T> = Result<T, Box<dyn Error + 'static>>;

/// A wrapper for any type implementing `Read` that counts the number of bytes read from the underlying buffer
/// on the inner reader.
///
/// The `ByteCounter` struct is used to track the number of bytes read from the inner reader, which is useful
/// when reading lines using the `.lines()` method. It provides a way to keep count of the total
/// bytes read and can be used to monitor the progress of reading operations.
///
/// Note the number of bytes read is Cumulative, and as such won't be reduced if the cursor is rewound.
///
/// # Example
///
/// ```rust,ignore
/// use std::io::{BufRead, Cursor};
/// use readfish_io::ByteCounter;
///
/// let data = "Line 1\nLine 2\nLine 3\n";
/// let cursor = Cursor::new(data);
/// let mut byte_counter = ByteCounter::new(cursor);
///
/// // Reading lines using the `.lines()` method will automatically update the byte count.
/// for line in byte_counter.lines() {
///     let line = line.unwrap();
///     println!("Line: {}", line);
/// }
///
/// println!("Total bytes read: {}", byte_counter.bytes_read());
/// ```
pub struct ByteCounter<R: Read> {
    /// The underlying Reader
    inner: R,
    /// Total number of bytes read.
    bytes_read: usize,
}

impl<R: Read> ByteCounter<R> {
    /// Creates a new [`ByteCounter`] instance wrapping the provided inner reader.
    ///
    /// This function constructs a new `ByteCounter` by taking ownership of the provided inner reader.
    /// The `ByteCounter` will keep track of the number of bytes read from the inner reader when
    /// consuming data.
    ///
    /// # Parameters
    ///
    /// - `inner`: The inner reader implementing the `Read` trait.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::io::Cursor;
    /// use readfish_io::ByteCounter;
    ///
    /// let data = b"Hello, world!";
    /// let cursor = Cursor::new(data);
    /// let byte_counter = ByteCounter::new(cursor);
    /// ```
    pub fn new(inner: R) -> Self {
        ByteCounter {
            inner,
            bytes_read: 0,
        }
    }
    /// Gets the total number of bytes read from the inner reader.
    ///
    /// This function returns the total number of bytes read from the inner reader since the creation
    /// of the `ByteCounter` instance. The count is updated every time data is read from the inner reader.
    ///
    /// # Returns
    ///
    /// The total number of bytes read so far.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::io::Cursor;
    /// use readfish_io::ByteCounter;
    ///
    /// let data = b"Hello, world!";
    /// let cursor = Cursor::new(data);
    /// let byte_counter = ByteCounter::new(cursor);
    ///
    /// // Reading from the inner reader will update the byte count
    /// let mut buffer = [0; 5];
    /// byte_counter.read_exact(&mut buffer).unwrap();
    ///
    /// let total_bytes_read = byte_counter.bytes_read();
    /// assert_eq!(total_bytes_read, 5);
    /// ```
    pub fn bytes_read(&self) -> usize {
        self.bytes_read
    }
}

impl<R: Read> Read for ByteCounter<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let bytes_read = self.inner.read(buf)?;
        self.bytes_read += bytes_read;
        Ok(bytes_read)
    }
}
impl<R: BufRead> BufRead for ByteCounter<R> {
    fn fill_buf(&mut self) -> ioResult<&[u8]> {
        // Delegate the call to the inner reader
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        // Keep track of the number of bytes read
        self.bytes_read += amt;
        // Delegate the call to the inner reader
        self.inner.consume(amt);
    }
}
// Create a new struct that wraps the GzDecoder and implements BufRead

/// A wrapper struct for a Gzip-compressed reader that also supports seeking.
/// It allows for consuming exact bytes to efficiently skip to a specific point in the reader.
/// I couldn't find the functionality I needed in the flate2 crate, so I had to implement it myself.
struct GzSeekable<R: Read + BufRead + Seek> {
    /// The inner reader, wrapping the GzDecoder
    inner: BufReader<GzDecoder<R>>,
}
impl<R: Read + BufRead + Seek> GzSeekable<R> {
    /// Create a new `GzSeekable` instance wrapping the provided reader `inner`.
    ///
    /// # Arguments
    ///
    /// * `inner`: The underlying reader implementing `Read`, `BufRead`, and `Seek`.
    ///
    /// # Returns
    ///
    /// A new `GzSeekable` instance.
    fn new(inner: R) -> Self {
        GzSeekable {
            inner: BufReader::new(GzDecoder::new(inner)),
        }
    }
}

impl<R: Read + BufRead + Seek> GzSeekable<R> {
    /// Consumes exactly `amt` bytes from the underlying Gzip-compressed reader.
    ///
    /// This method is used to efficiently skip to a specific point in the compressed data without
    /// decoding unnecessary bytes. It reads the data in blocks of 10,000 bytes (unless `amt` is smaller),
    /// and then reads the remaining bytes if there are any.
    ///
    /// # Arguments
    ///
    /// * `amt`: The number of bytes to consume from the reader.
    ///
    /// # Panics
    ///
    /// This function will panic if any I/O error occurs while reading from the underlying reader,
    /// which can happen if the input data is truncated or corrupted.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::io::{Read, Seek};
    /// use readfish_io::GzSeekable;
    /// use flate2::read::GzDecoder;
    ///
    /// // Create a Gzip-compressed reader from a Vec<u8>.
    /// let data = vec![31, 139, 8, 0, 0, 0, 0, 0, 2, 3, 103, 111, 95, 46, 98, 103, 7, 0, 40, 32, 13, 210, 73, 36, 0, 0, 0];
    /// let mut reader = GzSeekable::new(&data[..]);
    ///
    /// // Consume the first 10 bytes without decoding the remaining data.
    /// reader.consume_exact(10);
    /// ```
    ///
    /// In this example, the `consume_exact` method is used to skip the first 10 bytes of the compressed data
    /// without decoding the rest of the data, which can be useful for seeking to a specific position efficiently.

    fn consume_exact(&mut self, amt: usize) {
        let q = amt.checked_div_euclid(10000).unwrap();
        let r = amt.checked_rem_euclid(10000).unwrap();
        let mut buf = [0; 10000];
        for _i in 0..q {
            self.inner.read_exact(&mut buf).unwrap();
        }
        let mut buf = vec![0; r];
        self.inner.read_exact(&mut buf).unwrap();
    }
}

impl<R: Read + BufRead + Seek> Read for GzSeekable<R> {
    /// Read data from the Gzip-compressed reader.
    ///
    /// # Arguments
    ///
    /// * `buf`: A buffer to store the read data.
    ///
    /// # Returns
    ///
    /// The number of bytes read and stored in the buffer.
    fn read(&mut self, buf: &mut [u8]) -> ioResult<usize> {
        self.inner.read(buf)
    }
}

/// Get a buffered input reader from stdin or a file
///
/// # Arguments
///
/// * `path`: Optional path to a file. If provided, the function will attempt to open the file and create a buffered reader from it. If set to `None`, the function will create a buffered reader from stdin.
/// * `seek_bytes`: Optional number of bytes to seek into the file before reading. If `Some`, the function will seek to the specified number of bytes before reading the file. If `None`, the function will read the file from the beginning.
/// # Returns
///
/// A dynamic result that represents either a buffered reader on success or an error wrapped in a `Box<dyn Error + 'static>` on failure.
///
/// # Examples
///
/// ```rust,ignore
/// let reader = _get_reader_from_path(Some(PathBuf::from("path/to/file")))?;
/// ```
fn _get_reader_from_path(
    path: Option<PathBuf>,
    seek_bytes: Option<usize>,
) -> DynResult<Box<dyn BufRead + Send + 'static>> {
    let reader: Box<dyn BufRead + Send + 'static> = match path {
        Some(path) => {
            // stdin
            if path.as_os_str() == "-" {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, stdin()))
            } else {
                // open file
                let mut buf = BufReader::with_capacity(BUFFER_SIZE, File::open(path)?);
                buf.seek(SeekFrom::Start(seek_bytes.unwrap_or(0) as u64))?;
                Box::new(buf)
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
/// * `seek_bytes`: Optional number of bytes to seek into the file before reading. If `Some`, the function will seek to the specified number of bytes before reading the file. If `None`, the function will read the file from the beginning.
///
/// # Returns
///
/// A boxed trait object implementing `BufRead`, which can be used to read the contents of the file.
/// Uses the presence of a `.gz` or `.bgz` extension to decide
pub fn reader(
    filename: impl AsRef<Path>,
    seek_bytes: Option<usize>,
) -> Box<dyn BufRead + Send + 'static> {
    let ext = filename.as_ref().extension();
    let path: PathBuf = filename.as_ref().to_path_buf();
    // Handle Gzipped files first, since need to use flate2::read::GzDecoder

    if ext == Some(OsStr::new("gz")) {
        let file = match File::open(&path) {
            Err(why) => panic!("couldn't open {}: {}", path.display(), why),
            Ok(file) => file,
        };
        // Wrap the file in a `BufReader`
        let buf_reader = BufReader::new(file);

        // Wrap the `BufReader` in a `GzSeekable`, which are implementations that hopefully let us seek a certain amount into the file
        let mut wrapper = GzSeekable::new(buf_reader);
        // Read bytes indiscriminately, should be faster than iterating over lines?
        wrapper.consume_exact(seek_bytes.unwrap_or(0));
        let x = BufReader::new(wrapper);
        Box::new(x)

    // } else if ext == Some(OsStr::new("bgz")) {
    //     Box::new(BufReader::new(BgzfSyncReader::new(
    //         get_input(Some(path)).expect("Error: cannot read input file."),
    //     )))
    } else {
        _get_reader_from_path(Some(path), seek_bytes).expect("Error: cannot read input file")
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
fn _get_writer_from_path(path: Option<PathBuf>) -> DynResult<Box<dyn Write + Send + 'static>> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::{Cursor, Read};
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
    fn test_consume_exact_small() {
        // Create a Gzip-compressed data of size 50 bytes.
        let data = b"This is a Gzip-compressed data with 45 bytes.";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data).unwrap();
        let compressed_data = encoder.finish().unwrap();
        // Create a GzSeekable from the compressed data.
        let buf_reader = BufReader::new(Cursor::new(compressed_data));
        let mut reader = GzSeekable::new(buf_reader);

        // Consume 20 bytes.
        reader.consume_exact(20);
        let mut bc = ByteCounter::new(reader);
        bc.read_to_end(&mut Vec::new()).unwrap();
        // Verify that the reader position has been updated correctly.
        assert_eq!(bc.bytes_read(), 25);
    }

    #[test]
    fn test_consume_exact_large() {
        // Create a large Gzip-compressed data of size 100,000 bytes.
        let data = vec![0; 100_000];
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(&data).unwrap();
        let compressed_data = encoder.finish().unwrap();
        // Create a GzSeekable from the compressed data.
        let buf_reader = BufReader::new(Cursor::new(compressed_data));
        let mut reader = GzSeekable::new(buf_reader);
        // Consume 20 bytes.
        reader.consume_exact(43219);
        let mut bc = ByteCounter::new(reader);
        bc.read_to_end(&mut Vec::new()).unwrap();
        // Verify that the reader position has been updated correctly.
        assert_eq!(bc.bytes_read(), 56781);
    }

    #[test]
    fn test_byte_counter_new() {
        let data = b"Hello, world!";
        let cursor = Cursor::new(data);
        let byte_counter = ByteCounter::new(cursor);

        assert_eq!(byte_counter.bytes_read(), 0);
    }

    #[test]
    fn test_byte_counter_read() {
        let data = b"Hello, world!";
        let cursor = Cursor::new(data);
        let mut byte_counter = ByteCounter::new(cursor);

        let mut buffer = [0; 5];
        byte_counter.read_exact(&mut buffer).unwrap();

        assert_eq!(byte_counter.bytes_read(), 5);
    }
    #[test]
    fn test_byte_counter_read_past_end() {
        let data = b"Hello, world!";
        let cursor = Cursor::new(data);
        let mut byte_counter = ByteCounter::new(cursor);

        // Attempt to read more data than available in the cursor
        let mut buffer = [0; 100];
        let result = byte_counter.read(&mut buffer);

        assert_eq!(result.unwrap(), 13); // Only 13 bytes remaining in the cursor
        assert_eq!(byte_counter.bytes_read(), 13); // Total bytes read should be limited by cursor size
    }

    #[test]
    fn test_byte_counter_read_exact_past_end() {
        let data = b"Hello, world!";
        let cursor = Cursor::new(data);
        let mut byte_counter = ByteCounter::new(cursor);

        // Attempt to read more data than available in the cursor
        let mut buffer = [0; 100];
        let result = byte_counter.read_exact(&mut buffer);

        assert!(result.is_err()); // Should return an error as there's not enough data to fill the buffer
        assert_eq!(byte_counter.bytes_read(), 13); // Total bytes read should be limited by cursor size
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_reader() {
        let n_lines = reader(get_test_file("test_hum_4000.paf"), None)
            .lines()
            .count();
        assert_eq!(n_lines, 4148usize);

        let n_lines_gz = reader(get_test_file("test_hum_4000.paf.gz"), None)
            .lines()
            .count();
        assert_eq!(n_lines_gz, 4148usize);
    }
    #[test]
    fn test_reader_seek() {
        let mut reader = reader(get_test_file("test_hum_4000.paf"), Some(36));
        let mut line = String::new();
        let _line_read = reader.read_line(&mut line).unwrap();
        let line = line.trim();
        assert_eq!(line.len(), 119);
        assert_eq!(line, "635\t36\t632\t+\tNC_000007.14\t159345973\t115197340\t115197933\t453\t597\t60\ttp:A:P\tcm:i:66\ts1:i:452\ts2:i:63\tdv:f:0.0219\trl:i:138");
    }

    #[test]
    #[cfg_attr(miri, ignore)]

    fn test_reader_seek_gz() {
        let mut reader = reader(get_test_file("test_hum_4000.paf.gz"), Some(36));
        let mut line = String::new();
        let _line_read = reader.read_line(&mut line).unwrap();
        let line = line.trim();
        assert_eq!(line.len(), 119);
        assert_eq!(line, "635\t36\t632\t+\tNC_000007.14\t159345973\t115197340\t115197933\t453\t597\t60\ttp:A:P\tcm:i:66\ts1:i:452\ts2:i:63\tdv:f:0.0219\trl:i:138");
    }
}
