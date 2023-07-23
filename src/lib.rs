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
//! readfish_io - Custom functions and wrappers related IO functionality.
//! sequencing_summary - Sequencing summary related functionality.
mod channels;
pub mod nanopore;
mod paf;
pub mod readfish;
mod readfish_io;
mod sequencing_summary;
use std::{
    collections::HashMap,
    fmt,
    ops::Deref,
    path::{Path, PathBuf},
};

use itertools::Itertools;
use nanopore::format_bases;
use num_format::{Locale, ToFormattedString};
use paf::PafRecord;
use prettytable::{color, row, Attr, Cell, Row, Table};
use pyo3::prelude::*;
use readfish_io::DynResult;

/// Represents a summary of a contig or sequence from a sequencing experiment.
/// It includes various metrics related to the contig's characteristics and read mapping.
#[derive(Debug)]
pub struct ContigSummary {
    /// The name or identifier of the contig.
    pub name: String,
    /// The length of the contig in base pairs.
    pub length: usize,
    /// The mean read length of the mapped reads associated with this contig.
    pub mean_read_length: usize,
    /// The mean read quality of the mapped reads associated with this contig.
    pub mean_read_quality: f64,
    /// Yield of mapped reads
    pub total_bases: usize,
    /// The N50 metric for the contig, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length.
    pub n50: usize,
    /// The count of reads that are mapped on the target region (on-target reads).
    pub on_target_read_count: usize,
    /// The count of reads that are mapped off the target region (off-target reads).
    pub off_target_read_count: usize,
    /// The mean read length of on-target reads.
    pub mean_read_length_on_target: usize,
    /// The mean read length of off-target reads.
    pub mean_read_length_off_target: usize,
    /// The total yield (base pairs) of on-target reads for this contig.
    pub yield_on_target: usize,
    /// The total yield (base pairs) of off-target reads for this contig.
    pub yield_off_target: usize,
}
impl ContigSummary {
    /// Create a new `ContigSummary` instance with default values for all fields except `name` and `length`.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the contig.
    /// * `length` - The length of the contig.
    pub fn new(name: String, length: usize) -> Self {
        ContigSummary {
            name,
            length,
            mean_read_length: 0,
            mean_read_quality: 0.0,
            total_bases: 0,
            n50: 0,
            on_target_read_count: 0,
            off_target_read_count: 0,
            mean_read_length_on_target: 0,
            mean_read_length_off_target: 0,
            yield_on_target: 0,
            yield_off_target: 0,
        }
    }
    /// Get the total number of reads on the contig.
    pub fn total_reads(&self) -> usize {
        self.on_target_read_count + self.off_target_read_count
    }

    /// Mean read length of all reads on the contig.
    pub fn mean_read_length(&self) -> usize {
        self.total_bases / self.total_reads()
    }
    /// On target mean read length of all reads on the contig.
    pub fn on_target_mean_read_length(&self) -> usize {
        self.on_target_read_count / self.yield_on_target
    }
    /// Off target mean read length of all reads on the contig.
    pub fn off_target_mean_read_length(&self) -> usize {
        self.off_target_read_count / self.yield_off_target
    }
}
#[derive(Debug)]
/// Represents a summary of sequencing data, including various metrics related to the output of the experiment.
pub struct ConditionSummary {
    /// The name or identifier of the sequencing data.
    pub name: String,
    /// The total number of reads in the sequencing data.
    pub total_reads: usize,
    /// The count of reads that are mapped off the target regions (off-target reads).
    pub off_target_read_count: usize,
    /// The count of reads that are mapped to the target regions (on-target reads).
    pub on_target_read_count: usize,
    /// The percentage of off-target reads in the sequencing data.
    pub off_target_percent: f64,
    /// The total yield (base pairs) of off-target reads in the sequencing data.
    pub off_target_yield: usize,
    /// The total yield (base pairs) of on-target reads in the sequencing data.
    pub on_target_yield: usize,
    /// The mean read length of off-target reads.
    pub off_target_mean_read_length: usize,
    /// The mean read length of on-target reads.
    pub on_target_mean_read_length: usize,
    /// The mean read quality of off-target reads.
    pub off_target_mean_read_quality: f64,
    /// The mean read quality of on-target reads.
    pub on_target_mean_read_quality: f64,
    /// The N50 metric for the entire dataset, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length.
    pub n50: usize,
    /// The N50 metric for on-target reads, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length for on-target reads.
    pub on_target_n50: usize,
    /// The N50 metric for off-target reads, representing the length at which the cumulative
    /// sum of contig lengths reaches half of the total assembly length for off-target reads.
    pub off_target_n50: usize,
    /// A vector of `ContigSummary` representing summaries of individual contigs or sequences
    /// in the sequencing data.
    pub contigs: HashMap<String, ContigSummary>,
}

impl fmt::Display for ConditionSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Condition Name: {}", self.name)?;
        writeln!(f, "Total Reads: {}", self.total_reads)?;
        writeln!(f, "Off-Target Read Count: {}", self.off_target_read_count)?;
        writeln!(f, "On-Target Read Count: {}", self.on_target_read_count)?;
        writeln!(f, "Off-Target Percent: {:.2}%", self.off_target_percent)?;
        writeln!(f, "Off-Target Yield: {}", self.off_target_yield)?;
        writeln!(f, "On-Target Yield: {}", self.on_target_yield)?;
        writeln!(
            f,
            "Off-Target Mean Read Length: {}",
            self.off_target_mean_read_length
        )?;
        writeln!(
            f,
            "On-Target Mean Read Length: {}",
            self.on_target_mean_read_length
        )?;
        // writeln!(
        //     f,
        //     "Off-Target Mean Read Quality: {:.2}",
        //     self.off_target_mean_read_quality
        // )?;
        // writeln!(
        //     f,
        //     "On-Target Mean Read Quality: {:.2}",
        //     self.on_target_mean_read_quality
        // )?;
        // writeln!(f, "N50: {}", self.n50)?;
        // writeln!(f, "On-Target N50: {}", self.on_target_n50)?;
        // writeln!(f, "Off-Target N50: {}", self.off_target_n50)?;

        writeln!(f, "Contigs:")?;
        for (contig_name, contig_summary) in &self.contigs {
            writeln!(f, "  Contig Name: {}", contig_name)?;
            writeln!(f, "  Length: {}", contig_summary.length)?;
            // Print other fields from ContigSummary here
            // For example:
            // writeln!(f, "  Contig Mean Read Length: {}", contig_summary.mean_read_length)?;
        }
        Ok(())
    }
}

impl ConditionSummary {
    /// Update the `ConditionSummary` with information from the provided `PafRecord`.
    ///
    /// This method updates the fields of the `ConditionSummary` based on the information
    /// from the given `PafRecord`. It increments the appropriate read counts (on-target
    /// or off-target), calculates the mean read lengths and read qualities, updates the
    /// total reads count, and calculates the off-target percentage.
    ///
    /// # Arguments
    ///
    /// * `paf` - The [`PafRecord`] containing the information about the alignment.
    /// * `on_target` - A boolean flag indicating whether the alignment is on-target or off-target.
    ///
    /// # Returns
    ///
    /// This function returns a [`DynResult`] (a dynamic result that can contain any error).
    /// If the operation is successful, the `DynResult` will hold an `Ok(())`. Otherwise, it
    /// will hold an `Err` containing a helpful error message.
    pub fn update(&mut self, paf: PafRecord, on_target: bool) -> DynResult<()> {
        // update the condition struct
        self.total_reads += 1;
        if on_target {
            self.on_target_read_count += 1;
            self.on_target_yield += paf.query_length;
            self.on_target_mean_read_length = self.on_target_yield / self.on_target_read_count;
            // self.on_target_mean_read_quality += paf.tlen as f64;
        } else {
            self.off_target_read_count += 1;
            self.off_target_yield += paf.query_length;
            self.off_target_mean_read_length = self.off_target_yield / self.off_target_read_count;
            // self.off_target_mean_read_quality += paf.tlen as f64;
        }
        self.off_target_percent =
            self.off_target_read_count as f64 / self.total_reads as f64 * 100.0;
        let contig = self.get_or_add_contig(&paf.target_name, paf.target_length);
        contig.total_bases += paf.query_length;
        if on_target {
            contig.on_target_read_count += 1;
            contig.mean_read_length_off_target += paf.target_length;
            contig.mean_read_length_on_target =
                contig.yield_on_target / contig.on_target_read_count;
            // self.on_target_mean_read_quality += paf.tlen as f64;
        } else {
            contig.off_target_read_count += 1;
            contig.yield_off_target += paf.target_length;
            contig.mean_read_length_off_target =
                contig.yield_off_target / contig.off_target_read_count;
            // self.off_target_mean_read_quality += paf.tlen as f64;
        }
        // contig.mean_read_quality = paf.tlen;
        // contig.n50 = paf.tlen;
        // contig.on_target_read_count = paf.tlen;
        // contig.off_target_read_count = paf.tlen;

        Ok(())
    }
    /// Create a new `Summary` instance with default values for all fields except `name`.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the summary.
    pub fn new(name: String) -> Self {
        ConditionSummary {
            name,
            total_reads: 0,
            off_target_read_count: 0,
            on_target_read_count: 0,
            off_target_percent: 0.0,
            off_target_yield: 0,
            on_target_yield: 0,
            off_target_mean_read_length: 0,
            on_target_mean_read_length: 0,
            off_target_mean_read_quality: 0.0,
            on_target_mean_read_quality: 0.0,
            n50: 0,
            on_target_n50: 0,
            off_target_n50: 0,
            contigs: HashMap::new(),
        }
    }

    /// Get the name or identifier of the sequencing data.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Set the name or identifier of the sequencing data.
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    /// Get the total number of reads in the sequencing data.
    pub fn total_reads(&self) -> usize {
        self.total_reads
    }

    /// Set the total number of reads in the sequencing data.
    pub fn add_total_reads(&mut self, total_reads: usize) {
        self.total_reads += total_reads;
    }

    /// Get the count of reads that are mapped off the target regions (off-target reads).
    pub fn off_target_read_count(&self) -> usize {
        self.off_target_read_count
    }

    /// Set the count of reads that are mapped off the target regions (off-target reads).
    pub fn set_off_target_read_count(&mut self, off_target_read_count: usize) {
        self.off_target_read_count = off_target_read_count;
    }

    /// Get the count of reads that are mapped to the target regions (on-target reads).
    pub fn on_target_read_count(&self) -> usize {
        self.on_target_read_count
    }

    /// Set the count of reads that are mapped to the target regions (on-target reads).
    pub fn set_on_target_read_count(&mut self, on_target_read_count: usize) {
        self.on_target_read_count = on_target_read_count;
    }

    /// Get the percentage of off-target reads in the sequencing data.
    pub fn off_target_percent(&self) -> f64 {
        self.off_target_percent
    }

    /// Set the percentage of off-target reads in the sequencing data.
    pub fn set_off_target_percent(&mut self, off_target_percent: f64) {
        self.off_target_percent = off_target_percent;
    }

    /// Get the total yield (base pairs) of off-target reads in the sequencing data.
    pub fn off_target_yield(&self) -> usize {
        self.off_target_yield
    }

    /// Set the total yield (base pairs) of off-target reads in the sequencing data.
    pub fn set_off_target_yield(&mut self, off_target_yield: usize) {
        self.off_target_yield = off_target_yield;
    }

    /// Get the total yield (base pairs) of on-target reads in the sequencing data.
    pub fn on_target_yield(&self) -> usize {
        self.on_target_yield
    }

    /// Set the total yield (base pairs) of on-target reads in the sequencing data.
    pub fn set_on_target_yield(&mut self, on_target_yield: usize) {
        self.on_target_yield = on_target_yield;
    }
    /// Get the mean read length of off-target reads.
    pub fn off_target_mean_read_length(&self) -> usize {
        self.off_target_mean_read_length
    }

    /// Set the mean read length of off-target reads.
    pub fn set_off_target_mean_read_length(&mut self, off_target_mean_read_length: usize) {
        self.off_target_mean_read_length = off_target_mean_read_length;
    }

    /// Get the mean read length of on-target reads.
    pub fn on_target_mean_read_length(&self) -> usize {
        self.on_target_mean_read_length
    }

    /// Set the mean read length of on-target reads.
    pub fn set_on_target_mean_read_length(&mut self, on_target_mean_read_length: usize) {
        self.on_target_mean_read_length = on_target_mean_read_length;
    }

    /// Get the mean read quality of off-target reads.
    pub fn off_target_mean_read_quality(&self) -> f64 {
        self.off_target_mean_read_quality
    }

    /// Set the mean read quality of off-target reads.
    pub fn set_off_target_mean_read_quality(&mut self, off_target_mean_read_quality: f64) {
        self.off_target_mean_read_quality = off_target_mean_read_quality;
    }

    /// Get the mean read quality of on-target reads.
    pub fn on_target_mean_read_quality(&self) -> f64 {
        self.on_target_mean_read_quality
    }

    /// Set the mean read quality of on-target reads.
    pub fn set_on_target_mean_read_quality(&mut self, on_target_mean_read_quality: f64) {
        self.on_target_mean_read_quality = on_target_mean_read_quality;
    }

    /// Get the N50 metric for the entire dataset.
    pub fn n50(&self) -> usize {
        self.n50
    }

    /// Set the N50 metric for the entire dataset.
    pub fn set_n50(&mut self, n50: usize) {
        self.n50 = n50;
    }

    /// Get the N50 metric for on-target reads.
    pub fn on_target_n50(&self) -> usize {
        self.on_target_n50
    }

    /// Set the N50 metric for on-target reads.
    pub fn set_on_target_n50(&mut self, on_target_n50: usize) {
        self.on_target_n50 = on_target_n50;
    }

    /// Get the N50 metric for off-target reads.
    pub fn off_target_n50(&self) -> usize {
        self.off_target_n50
    }

    /// Set the N50 metric for off-target reads.
    pub fn set_off_target_n50(&mut self, off_target_n50: usize) {
        self.off_target_n50 = off_target_n50;
    }

    /// Get a reference to the vector of `ContigSummary`.
    pub fn contigs(&self) -> &HashMap<String, ContigSummary> {
        &self.contigs
    }

    /// Get a mutable reference to the vector of `ContigSummary`.
    pub fn contigs_mut(&mut self) -> &mut HashMap<String, ContigSummary> {
        &mut self.contigs
    }
    /// Get the ContigSummary associated with the given contig name or
    ///  add a new ContigSummary with the specified name and length if it doesn't exist.
    ///
    /// # Arguments
    ///
    /// * `contig` - The name of the contig.
    /// * `length` - The length of the contig.
    ///
    /// # Returns
    ///
    /// A reference to the ContigSummary associated with the contig name.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::collections::HashMap;
    /// use your_crate::ContigSummary;
    ///
    /// let mut contig_map: HashMap<String, ContigSummary> = HashMap::new();
    ///
    /// // Get an existing contig or add a new one if it doesn't exist
    /// let contig_name = "chr1".to_string();
    /// let contig_length = 1000;
    /// let contig_summary = contig_map.get_or_add_contig(contig_name.clone(), contig_length);
    ///
    /// // Now you can modify the ContigSummary fields or access its properties
    /// println!("Contig name: {}", contig_summary.name);
    /// println!("Contig length: {}", contig_summary.length);
    /// ```
    pub fn get_or_add_contig(&mut self, contig: &str, length: usize) -> &mut ContigSummary {
        self.contigs
            .entry(contig.to_string())
            .or_insert(ContigSummary::new(contig.to_string(), length))
    }

    /// get the total yield
    pub fn total_yield(&self) -> usize {
        self.on_target_yield + self.off_target_yield
    }
}

/// A struct representing a summary of conditions.
///
/// The `Summary` struct contains a hashmap where each key represents the name of a condition, and the corresponding value is a `ConditionSummary` struct
/// containing the summary information for that condition.
///
/// # Fields
///
/// * `conditions` - A hashmap containing the summary information for each condition. The key is a string representing the name of the condition,
/// and the value is a `ConditionSummary` struct containing the summary information for that condition.
///
/// # Examples
///
/// ```rust, ignore
/// use std::collections::HashMap;
/// use readfish_tools::{Summary, ConditionSummary};
///
/// // Create a new Summary
/// let mut summary = Summary {
///     conditions: HashMap::new(),
/// };
///
/// // Add some condition summaries
/// summary.conditions.insert(
///     "ConditionA".to_string(),
///     ConditionSummary {
///         // ... fill in the details for ConditionA ...
///     }
/// );
///
/// summary.conditions.insert(
///     "ConditionB".to_string(),
///     ConditionSummary {
///         // ... fill in the details for ConditionB ...
///     }
/// );
///
/// // Access a specific condition summary
/// if let Some(condition_summary) = summary.conditions.get("ConditionA") {
///     println!("Summary for ConditionA: {:?}", condition_summary);
/// }
/// ```
#[derive(Debug)]
pub struct Summary {
    /// Conditions summary for a given region or barcode.
    pub conditions: HashMap<String, ConditionSummary>,
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Summary:")?;
        // Todo rewrite to use Macro!
        for (condition_name, condition_summary) in &self.conditions {
            let mut condition_table = Table::new();
            condition_table.add_row(row!(bFg->format!("Condition Name: {}", condition_name)));
            condition_table.add_row(Row::new(vec![
                Cell::new("Total reads")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("# Off-target \nreads")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("# On-target \nreads")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Total Yield")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Off Target\n Yield")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("On Target\n yield")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
            ]));
            condition_table.add_row(Row::new(vec![
                // total reads
                Cell::new(
                    &condition_summary
                        .total_reads
                        .to_formatted_string(&Locale::en),
                )
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::GREEN)),
                // off target reads
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary
                        .off_target_read_count
                        .to_formatted_string(&Locale::en),
                    condition_summary.off_target_percent
                ))
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::GREEN)),
                // on target reads
                Cell::new(&format!(
                    "{} ({:.2}%)",
                    condition_summary
                        .on_target_read_count
                        .to_formatted_string(&Locale::en),
                    100_f64 - condition_summary.off_target_percent
                ))
                .with_style(Attr::Bold)
                .with_style(Attr::ForegroundColor(color::GREEN)),
                // total yield
                Cell::new(&format_bases(condition_summary.total_yield()))
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                // on target yield
                Cell::new(&format_bases(condition_summary.off_target_yield))
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                // on target yield
                Cell::new(&format_bases(condition_summary.on_target_yield))
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                // Cell::new(
                //     &condition_summary
                //         .off_target_mean_read_length
                //         .to_formatted_string(&Locale::en),
                // )
                // .with_style(Attr::Bold)
                // .with_style(Attr::ForegroundColor(color::GREEN)),
                // Cell::new(
                //     &condition_summary
                //         .on_target_mean_read_length
                //         .to_formatted_string(&Locale::en),
                // )
                // .with_style(Attr::Bold)
                // .with_style(Attr::ForegroundColor(color::GREEN)),
            ]));

            // writeln!(
            //     f,
            //     "  Off-Target Mean Read Quality: {:.2}",
            //     condition_summary.off_target_mean_read_quality
            // )?;
            // writeln!(
            //     f,
            //     "  On-Target Mean Read Quality: {:.2}",
            //     condition_summary.on_target_mean_read_quality
            // )?;
            // writeln!(f, "  N50: {}", condition_summary.n50)?;
            // writeln!(f, "  On-Target N50: {}", condition_summary.on_target_n50)?;
            // writeln!(f, "  Off-Target N50: {}", condition_summary.off_target_n50)?;
            condition_table.printstd();
            writeln!(f, "Contigs:")?;
            let mut contig_table = Table::new();
            // Create a custom format with left-leading spaces
            contig_table.get_format().indent(10);
            contig_table.add_row(Row::new(vec![
                Cell::new("Contig")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Contig Length")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Read count")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Yield")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Mean \nRead Length")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("On Target\n Reads")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                Cell::new("Off \nTarget Reads")
                    .with_style(Attr::Bold)
                    .with_style(Attr::ForegroundColor(color::GREEN)),
            ]));
            for (contig_name, contig_summary) in condition_summary
                .contigs
                .iter()
                .sorted_by(|(key1, _), (key2, _)| natord::compare(key1, key2))
            {
                contig_table.add_row(Row::new(vec![
                    Cell::new(contig_name)
                        .with_style(Attr::Bold)
                        .with_style(Attr::ForegroundColor(color::GREEN)),
                    Cell::new(&contig_summary.length.to_formatted_string(&Locale::en))
                        .with_style(Attr::ForegroundColor(color::GREEN)),
                    Cell::new(
                        &contig_summary
                            .total_reads()
                            .to_formatted_string(&Locale::en),
                    )
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                    Cell::new(&contig_summary.total_bases.to_formatted_string(&Locale::en))
                        .with_style(Attr::ForegroundColor(color::GREEN)),
                    Cell::new(
                        &contig_summary
                            .mean_read_length()
                            .to_formatted_string(&Locale::en),
                    )
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                    Cell::new(
                        &contig_summary
                            .on_target_read_count
                            .to_formatted_string(&Locale::en),
                    )
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                    Cell::new(
                        &contig_summary
                            .off_target_read_count
                            .to_formatted_string(&Locale::en),
                    )
                    .with_style(Attr::ForegroundColor(color::GREEN)),
                ]));
                // Print other fields from ContigSummary here
                // For example:
                // writeln!(f, "    Contig Mean Read Length: {}", contig_summary.mean_read_length)?;
            }
            contig_table.printstd();
        }
        Ok(())
    }
}

impl Summary {
    /// Create a new `Summary` instance with default values for all fields.
    fn new() -> Self {
        Summary {
            conditions: HashMap::new(),
        }
    }

    /// Get the summary for the specified condition. If the condition does not exist in the
    /// `Summary`, it will be created with default values.
    ///
    /// # Arguments
    ///
    /// * `condition_name`: A `String` representing the name or identifier of the condition.
    ///
    /// # Returns
    ///
    /// A reference to the `ConditionSummary` for the specified condition.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use std::collections::HashMap;
    ///
    /// let mut summary = Summary::new();
    ///
    /// // Get or add the condition with the name "Condition A"
    /// let condition_a = summary.conditions("Condition A".to_string());
    ///
    /// // Modify the fields of the condition summary
    /// condition_a.set_total_reads(10000);
    /// condition_a.set_on_target_read_count(8000);
    /// condition_a.set_off_target_read_count(2000);
    /// // ...
    ///
    /// // Get or add another condition
    /// let condition_b = summary.conditions("Condition B".to_string());
    /// // ...
    /// ```
    pub fn conditions<T: Deref<Target = str>>(
        &mut self,
        condition_name: T,
    ) -> &mut ConditionSummary {
        self.conditions
            .entry(condition_name.to_string())
            .or_insert(ConditionSummary::new(condition_name.to_string()))
    }
}
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
pub fn _demultiplex_paf(
    toml_path: impl AsRef<Path>,
    paf_path: impl AsRef<Path>,
    sequencing_summary_path: Option<impl AsRef<Path>>,
    print_summary: bool,
    _csv_out: Option<impl AsRef<Path>>,
) {
    let toml_path = toml_path.as_ref();
    let paf_path = paf_path.as_ref();
    let toml = readfish::Conf::from_file(toml_path);
    let mut paf = paf::Paf::new(paf_path);
    let seq_sum =
        sequencing_summary_path.map(|path| sequencing_summary::SeqSum::from_file(path).unwrap());
    let mut seq_sum = seq_sum;
    let mut summary = Summary::new();
    paf.demultiplex(&toml, seq_sum.as_mut(), Some(&mut summary))
        .unwrap();
    if print_summary {
        println!("{}", summary);
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn demultiplex_paf(toml_path: PathBuf, paf_path: PathBuf, seq_sum_path: PathBuf) -> PyResult<()> {
    _demultiplex_paf(
        toml_path,
        paf_path,
        Some(seq_sum_path),
        true,
        None::<String>,
    );
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn readfish_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(demultiplex_paf, m)?)?;
    Ok(())
}
