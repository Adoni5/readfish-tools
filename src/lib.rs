#![deny(missing_docs)]
#![warn(clippy::missing_docs_in_private_items)]
#![allow(dead_code)]
//! # Readfish-tools
//!
//! `readfish-tools` is a collection of utilities to provide a standardised way of analysing
//! readfish runs that have been run. Currently the accepted analysable inputs are sequencing summary files,
//! BAM of all produced FASTQ, and the `TOML` file that was used to configure the readfish run.
//!
//! The intention is to demultiplex a bam/paf/sequencing summary into regions and barcodes then have methods to provide the
//! summary stats for this function.
//!
//! The crate is split into modules handling separate functionalities.
//!
//! ## Modules
//! nanopore - Flowcell related functionality.
//! channels - Channel Hashmaps for MinION and Flongle.
//! paf - PAF related functionality.
//! readfish - Readfish TOML related functionality.
//!
//!
mod channels;
pub mod nanopore;
mod paf;
pub mod readfish;
mod readfish_io;
mod sequencing_summary;
use std::path::Path;

/// Demultiplex PAF records based on the specified configuration.
///
/// This function takes two file paths as inputs, `toml_path` and `paf_path`, representing
/// the paths to the TOML configuration file and the PAF file, respectively. The TOML configuration
/// is read using the `readfish::Conf::from_file` function, and the PAF file is opened and checked using the
/// `paf::open_paf_for_reading` function. The resulting PAF records are then demultiplexed based on the
/// information provided in the configuration file.
///
/// Note: The current implementation initializes a new `paf::Paf` object with a hardcoded PAF file
/// path "resources/test_paf_With_seq_sum.paf" and calls its `demultiplex` method with the parsed
/// TOML configuration. However, the line is commented out, so the actual demultiplexing process
/// is not performed. Please ensure that the proper PAF object is used and uncommented to perform
/// the demultiplexing.
///
/// If there are barcodes present in the Conf TOML file, and the barcode_arrangement column is missing from the
/// the sequencing summary file, the function will panic.
///
/// # Arguments
///
/// * `toml_path`: The file path to the TOML configuration file.
/// * `paf_path`: The file path to the PAF file to be demultiplexed.
///
/// # Examples
///
/// ```rust,ignore
/// use std::path::Path;
/// demultiplex_paf("config.toml", "file.paf");
/// ```
///
pub fn demultiplex_paf(
    toml_path: impl AsRef<Path>,
    paf_path: impl AsRef<Path>,
    sequencing_summary_path: Option<String>,
) {
    let toml_path = toml_path.as_ref();
    let paf_path = paf_path.as_ref();
    let toml = readfish::Conf::from_file(toml_path);
    let mut paf = paf::Paf::new(paf_path);
    let seq_sum =
        sequencing_summary_path.map(|path| sequencing_summary::SeqSum::from_file(path).unwrap());
    let mut seq_sum = seq_sum;
    paf.demultiplex(&toml, seq_sum.as_mut()).unwrap();
}
