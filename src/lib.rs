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

use csv::ReaderBuilder;
use serde::Deserialize;
use std::{
    collections::HashMap,
    hash::{Hash, Hasher},
    io::Cursor,
    path::{Path, PathBuf},
};
use toml::{map::Map, Table, Value};

mod channels;
pub mod nanopore;

/// Action types that can be taken once a decision (one of single_off, single_on, multi_off, multi_on, no_map, no_seq, exceeded_max_chunks, below_min_chunks)
/// has been made.
#[derive(Debug)]
enum Action {
    /// Read would be unblocked
    Unblock,
    /// Complete sequencing naturally
    StopReceiving,
    /// Proceed with sequencing
    Proceed,
}
/// Type for the Contig -> coordinates hashmap.
type HashedTargets = HashMap<String, Vec<(usize, usize)>>;

impl From<&str> for Action {
    fn from(source: &str) -> Action {
        match source {
            "unblock" => Action::Unblock,
            "stop_receiving" => Action::StopReceiving,
            "proceed" => Action::Proceed,
            _ => {
                panic!("Unknown Action given")
            }
        }
    }
}

/// The _Condition struct holds the settings lifted from the TOML file, for each
/// region of the flowcell or barcode.
#[derive(Debug)]
struct _Condition {
    /// The name of the Condition (Barcode/Region).
    name: String,
    /// The minimum number of read chunks that have to be captured for a read to be processed. Default if not met is to proceed.
    min_chunks: u8,
    /// The maximum number of read chunks that can be captured for a read. Default if exceed is to unblock.
    max_chunks: u8,
    /// The targets associated with the Condition.
    targets: Targets,
    /// The action to perform when an alignment returns one single primary mapping, outside of any target regions.
    single_off: Action,
    /// The action to perform when an alignment returns one single primary mapping, inside of a target regions.
    single_on: Action,
    /// The action to perform when an alignment returns multiple primary mappings, all outside of any target regions.
    multi_off: Action,
    /// The action to perform when an alignment returns multiple primary mappings, at LEAST ONE of which is inside of a target region.
    multi_on: Action,
    /// The action to perform when no alignments are returned for this read.
    no_map: Action,
    /// The action to perform when no sequence is produced for this read sequence.
    no_seq: Action,
}

#[derive(Debug)]
/// Represents a region of the flow cell, denoted in the configuration toml as
///
/// ```toml
///
///    [[regions]]
///    name = "Rapid_CNS"
///    min_chunks = 1
///    max_chunks = 4
///    targets = "resources/panel_adaptive_nogenenames_20122021_hg38.bed"
///    single_off = "unblock"
///    multi_off = "unblock"
///    single_on = "stop_receiving"
///    multi_on = "stop_receiving"
///    no_seq = "proceed"
///    no_map = "proceed"
/// ```
/// All the parsed fields are stored with a _Condition struct, as they could also be from a barcodes table.
struct Region {
    /// The parsed region settings.
    condition: _Condition,
}

/// Represents a barcode on the sequencing library. This supercedes any regions.
///
/// ```toml
///
//[barcodes.barcode02]
//name = "barcode02"
//control = false
//min_chunks = 0
//max_chunks = 4
//targets = []
//single_on = "unblock"
//multi_on = "unblock"
//single_off = "unblock"
//multi_off = "unblock"
//no_seq = "proceed"
///no_map = "unblock"
/// ```
///
/// All the parsed fields are stored with a _Condition struct, as they could also be from a regions table.
#[derive(Debug)]
struct Barcode {
    /// The parsed barcode settings.
    condition: _Condition,
}

impl From<&Map<String, Value>> for _Condition {
    fn from(source: &Map<String, Value>) -> Self {
        let targets: TargetType = source.get("targets").unwrap().into();
        let target: Targets = Targets::new(targets);
        _Condition {
            name: source.get("name").unwrap().as_str().unwrap().to_string(),
            min_chunks: source
                .get("min_chunks")
                .unwrap_or(&toml::Value::Integer(0))
                .as_integer()
                .unwrap()
                .try_into()
                .unwrap(),

            max_chunks: source
                .get("min_chunks")
                .unwrap_or(&toml::Value::Integer(4))
                .as_integer()
                .unwrap()
                .try_into()
                .unwrap(),
            targets: target,
            single_off: source.get("single_off").unwrap().as_str().unwrap().into(),
            single_on: source.get("single_on").unwrap().as_str().unwrap().into(),
            multi_on: source.get("multi_on").unwrap().as_str().unwrap().into(),
            multi_off: source.get("multi_off").unwrap().as_str().unwrap().into(),
            no_map: source.get("no_map").unwrap().as_str().unwrap().into(),
            no_seq: source.get("no_seq").unwrap().as_str().unwrap().into(),
        }
    }
}

/// Strand that the target is on.
#[derive(Debug, Hash, PartialEq)]
enum Strand {
    /// Represents he forward (sense) strand
    Forward,
    /// Represents he reverse (anti-sense) strand
    Reverse,
}

/// A wrapper for the Strand, which implements Hash and Eq, allowing the Strand enum to be used for
/// a HashMap key.
///
/// Implements to_string and AsRef str to get string representations, so we can take it along with multiple other types into functions
/// that need the strand.
#[derive(PartialEq, Debug)]
struct StrandWrapper(Strand);

impl Eq for StrandWrapper {}

impl Hash for StrandWrapper {
    fn hash<H: Hasher>(&self, state: &mut H) {
        std::mem::discriminant(&self.0).hash(state);
    }
}

impl From<&str> for Strand {
    fn from(source: &str) -> Strand {
        match source {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            "1" => Strand::Forward,
            "-1" => Strand::Reverse,
            _ => Strand::Forward,
        }
    }
}

impl ToString for Strand {
    fn to_string(&self) -> String {
        match self {
            Strand::Forward => "+".to_string(),
            Strand::Reverse => "-".to_string(),
        }
    }
}

impl AsRef<str> for Strand {
    fn as_ref(&self) -> &str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        }
    }
}
/// TargetRype Enum, represents whther targets were listed directly in the TOML file
/// or a path to a targets containing file was given.
#[derive(Clone, Debug)]
enum TargetType {
    /// Variant representing targets that were given directly in the TOML file.
    Direct(Vec<String>),
    /// Variant representing targets that were given as a path to a file that contains targets.
    ViaFile(PathBuf),
}
/// Represents a BED record, which is read from a BedFILE. All six columns are expected, however we do not use _name or _score.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct BedRecord {
    /// The contig or chromosome name associated with the record.
    contig: String,
    /// The starting position of the record.
    start: usize,
    /// The stopping position of the record.
    stop: usize,
    /// The name associated with the record (optional).
    _name: String,
    /// The score associated with the record (optional).
    _score: String,
    /// The strand information of the record.
    strand: String,
}

/// CSV record parsed from targets specified in TOML file,
/// If A bed file is provided, the six records are taken and placed in a
/// BedRecord. This BedRecord is then converted into a CsvRecord.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CsvRecord {
    /// Contig the target is on
    pub contig: String,
    /// Optional start coordinate of target
    #[serde(default)]
    pub start: Option<usize>,
    /// Optional stop coordinate of target. Required if start is present
    #[serde(default)]
    pub stop: Option<usize>,
    /// Optional strand target is on. .One of "+"/"-". Required if Start/Stop are provided.
    #[serde(default)]
    pub strand: Option<String>,
}

impl From<BedRecord> for CsvRecord {
    fn from(source: BedRecord) -> CsvRecord {
        CsvRecord {
            contig: source.contig,
            start: Some(source.start),
            stop: Some(source.stop),
            strand: Some(source.strand),
        }
    }
}

impl CsvRecord {
    /// Checks if the structure has valid coordinates.
    ///
    /// Returns `true` if both `start` and `stop` fields have values,
    /// indicating that the structure has valid coordinates. Otherwise, returns `false`.
    ///
    /// # Examples
    ///
    /// ```rust, ignore
    /// # use readfish_tools::CsvRecord;
    ///
    /// let record = CsvRecord {
    ///     contig: "chr1".to_string(),
    ///     start: Some(100),
    ///     stop: Some(200),
    ///     strand: Some("+".to_string()),
    /// };
    ///
    /// assert!(record.has_coords());  // Returns true
    ///
    /// let invalid_record = CsvRecord {
    ///     contig: "chr2".to_string(),
    ///     start: None,
    ///     stop: None,
    ///     strand: Some("-".to_string()),
    /// };
    ///
    /// assert!(!invalid_record.has_coords());  // Returns false
    /// ```
    fn has_coords(&self) -> bool {
        self.start.is_some() && self.stop.is_some()
    }

    /// Retrieves the coordinates from the structure.
    ///
    /// Returns a tuple containing the start and stop coordinates of the structure.
    /// If the structure has valid coordinates (i.e., `has_coords()` is true),
    /// the actual values of the `start` and `stop` fields are returned.
    /// Otherwise, a default range of (0, usize::MAX) is returned.
    ///
    /// # Examples
    ///
    /// ```rust, ignore
    /// # use readfish_tools::CsvRecord;
    ///
    /// let record = CsvRecord {
    ///     contig: "chr1".to_string(),
    ///     start: Some(100),
    ///     stop: Some(200),
    ///     strand: Some("+".to_string()),
    /// };
    ///
    /// assert_eq!(record.get_coords(), (100, 200));
    ///
    /// let invalid_record = CsvRecord {
    ///     contig: "chr2".to_string(),
    ///     start: None,
    ///     stop: Some(300),
    ///     strand: Some("-".to_string()),
    /// };
    ///
    /// assert_eq!(invalid_record.get_coords(), (0, usize::MAX));
    /// ```
    fn get_coords(&self) -> (usize, usize) {
        if self.has_coords() {
            (self.start.unwrap(), self.stop.unwrap())
        } else {
            (0, usize::MAX)
        }
    }

    /// Retrieves the strand information from the struct.
    ///
    /// This function returns an `Option<Strand>` representing the strand information stored in the struct.
    /// If the `strand` field is [`Some`], the function maps the string value to the corresponding [`Strand`] enum variant
    /// using the `from` method. If the `strand` field is [`None`], the function returns `None`.
    ///
    /// # Returns
    ///
    /// An `Option<Strand>` representing the strand information, or `None` if no strand is available.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// # use readfish_tools::CsvRecord;
    ///
    /// let record = CsvRecord {
    ///     contig: "chr1".to_string(),
    ///     start: Some(100),
    ///     stop: Some(200),
    ///     strand: Some("+".to_string()),
    /// };
    /// let strand = record.get_strand();
    /// assert_eq(strand, Some(Strand::Forward))
    /// ```
    fn get_strand(&self) -> Option<Strand> {
        self.strand
            .as_ref()
            .map(|strand_string| Strand::from(strand_string.as_str()))
    }
}

impl From<&Value> for TargetType {
    fn from(source: &Value) -> TargetType {
        match source.is_array() {
            true => TargetType::Direct(
                source
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_str().unwrap().to_string())
                    .collect(),
            ),
            false => TargetType::ViaFile(PathBuf::from(source.as_str().unwrap())),
        }
    }
}

/// Represents a configuration for a flowcell.
#[derive(Debug)]
struct Conf {
    /// The total number of channels on the flowcell.
    channels: usize,
    /// The regions of the flowcell. contains the name of the region and the Action to take for each Alignment type.
    regions: Vec<Region>,
    /// The barcodes from the sequencing library.
    barcodes: Vec<Barcode>,
    /// The mapping of channel number to the index of the region that channel belongs to.
    _channel_map: HashMap<u16, u8>,
}
#[derive(Debug)]
/// Holds the targets for a given region or barcode.
struct Targets {
    /// The target string as listed in the Toml. Can either be an array of strings, in which case that is assumed to be the targets themselves, or a string,
    /// which is assumed to be a file path to a file containing the targets.
    value: TargetType,
    /// A hashmap containg the targets themselves, in the form of
    /// Strand => Contig => Start and stop target coordinates.
    _targets: HashMap<StrandWrapper, HashedTargets>,
}

impl Targets {
    /// Creates a new instance of [`Targets`] with the provided target data.
    ///
    /// This function takes the target data in the form of [`TargetType`] and constructs a new [`Targets`] struct.
    /// The [`TargetType`] can either be an array of strings, representing the targets themselves,
    /// or a string representing a file path to a file containing the targets (.bed or .csv).
    ///
    /// If the targets is an array of strings, they must be im the format "contig, start, stop, strand" OR "contig".
    /// For example:
    ///
    /// ```toml
    /// targets = ["chr2,10.20.+", "chr1"]
    /// ```
    ///
    /// If only the contig is provided, it is assumed that the whole contig is the target, on BOTH strands. See example below.
    ///
    /// The target data is stored in the `value` field, while the parsed targets are stored in the `_targets` field
    /// as a hashmap with strand, contig, and start/stop target coordinates.
    ///
    /// # Arguments
    ///
    /// * `targets` - The target data in the form of [`TargetType`].
    ///
    /// # Examples
    ///
    /// ```rust, ignore
    /// # use my_module::{Targets, TargetType};
    ///
    /// let target_data = TargetType::Direct(vec!["chr1".to_string(), "chr2,10.20.+".to_string()]);
    /// let targets = Targets::new(target_data);
    ///
    /// assert_eq!(targets.value, TargetType::Array(vec!["chr1".to_string(), "chr2,10.20.+".to_string()]));
    /// assert_eq!(targets._targets.len(), 2);
    ///
    /// println!("{:#?}", targets._targets)
    /// // {
    /// //    StrandWrapper(Forward): {"chr1": [(0, 18_446_744_073_709_551_615)], "chr2": [(10, 20)]}
    /// //    StrandWrapper(Reverse): {"chr1": [(0, 18_446_744_073_709_551_615)]}
    /// // }
    /// // NOTE the single contig target chr1 is on both strands in its entirety.
    ///
    /// ```
    fn new(targets: TargetType) -> Targets {
        let t = targets.clone();
        Targets {
            value: targets,
            _targets: Targets::from_parsed_toml(t),
        }
    }

    /// Inserts target coordinates into the `targets` hashmap based on the provided record and strand.
    ///
    /// This function takes a mutable reference to the `targets` hashmap, a reference to a [`CsvRecord`],
    /// and a variant from the `strand` Enum. It inserts the record coordinates into a Vec at the lowest level of
    /// the `targets` hashmap based on the strand and contig.
    ///
    /// If the strand does not exist in the `targets` hashmap, a new entry is created for the strand,
    /// and an empty hashmap is inserted for the contig. If the contig does not exist for the strand,
    /// a new entry is created for the contig, and an empty vector is inserted to store the coordinates.
    ///
    /// The record coordinates are retrieved using the `get_coords()` method from the [`CsvRecord`] struct.
    ///
    /// # Arguments
    ///
    /// * `targets` - A mutable reference to the `HashMap<StrandWrapper, HashedTargets>]` where the record will be inserted.
    /// * `record` - A reference to the `CsvRecord` containing the record information.
    /// * `strand` - The strand information associated with the record.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use readfish_tools::{insert_into_targets, CsvRecord, Strand, StrandWrapper, HashedTargets};
    /// use std::collections::HashMap;
    ///
    /// let mut targets: HashMap<StrandWrapper, HashedTargets> = HashMap::new();
    ///
    /// let record = CsvRecord {
    ///     contig: "chr1".to_string(),
    ///     start: Some(100),
    ///     stop: Some(200),
    ///     strand: Some("+".to_string()),
    /// };
    ///
    /// insert_into_targets(&mut targets, &record, "+");
    ///
    /// assert_eq!(targets.len(), 1);
    /// assert_eq!(targets.get(&StrandWrapper(Strand::Forward)).unwrap().len(), 1);
    /// assert_eq!(targets.get(&StrandWrapper(Strand::Forward)).unwrap().get("chr1").unwrap().len(), 1);
    /// assert_eq!(targets.get(&StrandWrapper(Strand::Forward)).unwrap().get("chr1").unwrap()[0], (100, 200));
    /// ```
    fn insert_into_targets(
        targets: &mut HashMap<StrandWrapper, HashedTargets>,
        record: &CsvRecord,
        strand: Strand,
    ) {
        let coords = targets
            .entry(StrandWrapper(strand))
            .or_insert(HashMap::new())
            .entry(record.contig.clone())
            .or_insert(Vec::with_capacity(1000));
        coords.push(record.get_coords())
    }

    /// Creates a hashmap of targets from the parsed TOML data.
    ///
    /// This function takes the `targets` data in the form of [`TargetType`]] and constructs a hashmap of targets
    /// grouped by strand and contig, with start and stop coordinates as values. The `targets` can be provided
    /// either as a direct array of target strings or as a path to a CSV or BED file containing the targets.
    ///
    /// If `targets` is of type [`TargetType::Direct`], the function treats the data as direct target strings,
    /// parses them as CSV data, and populates the hashmap with the targets grouped by strand and contig.
    /// If `targets` is of type [`TargetType::ViaFile`], the function treats the data as a file path,
    /// determines the file type (CSV or BED), and parses the data accordingly to populate the hashmap.
    ///
    /// The function uses the [`CsvRecord`] struct for deserialization of CSV records and the [`BedRecord`] struct
    /// for deserialization of BED records. The appropriate deserialization is performed based on the file type.
    ///
    /// After populating the hashmap, the function merges overlapping intervals within each contig
    /// using the [`Self::_merge_intervals()`] helper function.
    ///
    /// # Arguments
    ///
    /// * `targets` - The target data in the form of [`TargetType`].
    ///
    /// # Returns
    ///
    /// A hashmap of targets grouped by strand and contig, with start and stop coordinates as values.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use readfish_tools::{Targets::from_parsed_toml, TargetType, CsvRecord};
    /// use std::collections::HashMap;
    /// use std::path::PathBuf;
    ///
    /// let target_data = TargetType::Direct(vec![
    ///     "chr2,".to_string(),
    ///     "chr1,10,20,+".to_string(),
    /// ]);
    ///
    /// let targets = from_parsed_toml(target_data);
    ///
    /// assert_eq!(targets.len(), 2);
    /// assert_eq!(targets.get(&StrandWrapper(Strand::Forward)).unwrap().get("chr2").unwrap()[0], (0_usize, usize::MAX));
    /// assert_eq!(targets.get(&StrandWrapper(Strand::Forward)).unwrap().get("chr1").unwrap()[0], (10_usize,20_usize));
    /// ```
    fn from_parsed_toml(
        targets: TargetType,
    ) -> HashMap<StrandWrapper, HashMap<String, Vec<(usize, usize)>>> {
        let mut results = HashMap::new();
        let mut bed_file = false;
        let mut delim = b',';
        match targets {
            TargetType::Direct(target_vec) => {
                if target_vec.is_empty() {
                    return results;
                }
                let csv_data = target_vec.join("\n");
                let file = Cursor::new(csv_data);
                let mut reader = ReaderBuilder::new()
                    .flexible(true)
                    .has_headers(false)
                    .delimiter(delim)
                    .from_reader(file);
                for record in reader.records() {
                    let record = record.unwrap();
                    let record: CsvRecord = record.deserialize(None).unwrap();
                    if record.has_coords() {
                        Targets::insert_into_targets(
                            &mut results,
                            &record,
                            record.get_strand().unwrap(),
                        );
                    } else {
                        Targets::insert_into_targets(&mut results, &record, Strand::Forward);
                        Targets::insert_into_targets(&mut results, &record, Strand::Reverse);
                    }
                }
            }
            TargetType::ViaFile(file_path) => {
                // TODO won't handle gzipped bed files
                if file_path.extension().unwrap() == "bed" {
                    bed_file = true;
                    delim = b'\t';
                }
                let mut rdr = ReaderBuilder::new()
                    .delimiter(delim)
                    .flexible(true)
                    .has_headers(false)
                    .from_path(file_path)
                    .unwrap();
                for record in rdr.records() {
                    let record = record.unwrap();
                    let record: CsvRecord = match bed_file {
                        true => {
                            let x: BedRecord = record.deserialize(None).unwrap();
                            x.into()
                        }
                        false => {
                            let x: CsvRecord = record.deserialize(None).unwrap();
                            x
                        }
                    };
                    // Has coordinates and strand provided
                    if record.has_coords() {
                        Targets::insert_into_targets(
                            &mut results,
                            &record,
                            record.get_strand().unwrap(),
                        );
                    }
                }
            }
        }
        results.iter_mut().for_each(|(_strand, contig_hashmap)| {
            contig_hashmap
                .iter_mut()
                .for_each(|(_, v)| *v = Targets::_merge_intervals(v))
        });
        results
    }

    /// Merges overlapping intervals within a vector of intervals.
    ///
    /// This function takes a mutable reference to a vector of intervals represented as tuples `(usize, usize)`
    /// and merges any overlapping intervals into collapsed ranges. The intervals are expected to be sorted
    /// based on the starting index before calling this function.
    ///
    /// If the number of intervals is less than 2, the function returns a clone of the input vector as there
    /// are no overlapping intervals to merge.
    ///
    /// The function iterates over the sorted intervals and maintains a current range. For each interval,
    /// if it overlaps with the current range, the end index of the current range is updated to the maximum
    /// of the current end index and the interval's end index. If the interval is non-overlapping, the
    /// current range is added to the collapsed ranges and updated to the new interval. If it's the first
    /// range encountered, the current range is initialized. Finally, the last current range (if any) is added
    /// to the collapsed ranges.
    ///
    /// The resulting collapsed ranges are returned as a new vector.
    ///
    /// # Arguments
    ///
    /// * `intervals` - A mutable reference to a vector of intervals to be merged.
    ///
    /// # Returns
    ///
    /// A vector of collapsed ranges after merging overlapping intervals.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    ///
    /// let mut intervals = vec![(1, 5), (4, 9), (10, 15), (13, 18)];
    /// let collapsed_ranges = Targets::_merge_intervals(&mut intervals);
    ///
    /// assert_eq!(collapsed_ranges, vec![(1, 9), (10, 18)]);
    /// ```
    fn _merge_intervals(intervals: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        // ToDo consider explicitly forbidding start > end or end < start
        let n_args = intervals.len();
        if n_args < 2 {
            return intervals.clone();
        }
        intervals.sort(); // Sort the ranges based on the starting index
        let mut collapsed_ranges: Vec<(usize, usize)> = Vec::new();
        let mut current_range: Option<(usize, usize)> = None;
        for &(start, end) in intervals.iter() {
            if let Some((current_start, current_end)) = current_range {
                if start <= current_end {
                    // Overlapping range, update the current range's end index
                    current_range = Some((current_start, current_end.max(end)));
                } else {
                    // Non-overlapping range, add the current range and update the current range
                    collapsed_ranges.push((current_start, current_end));
                    current_range = Some((start, end));
                }
            } else {
                // First range encountered, initialize the current range
                current_range = Some((start, end));
            }
        }
        // Add the last current range (if any)
        if let Some((current_start, current_end)) = current_range {
            collapsed_ranges.push((current_start, current_end));
        }
        collapsed_ranges
    }

    /// Checks if the given coordinate falls within any of the target intervals for the specified contig and strand.
    ///
    /// This function takes a reference to a [`CsvRecord`] struct and performs a lookup in the [`Targets`] struct's
    /// `_targets` hashmap to retrieve the intervals for the specified contig and strand. It then checks if the
    /// given coordinate falls within any of the target intervals by iterating over the intervals and performing
    /// the comparison.
    ///
    /// The function expects the `strand` argument to implement the [`ToString`] trait, which allows the function
    /// to convert it to a [`String`]. The `strand` is then converted to the [`Strand`] enum type using the `into()`
    /// method.
    ///
    /// # Generic Parameters
    ///
    /// * `T` - The type of the `strand` argument that implements the [`ToString`] trait.
    ///
    /// # Arguments
    ///
    /// * `contig` - The contig string to lookup the intervals for.
    /// * `strand` - The strand value to lookup the intervals for. It is expected to be convertible to a [`String`].
    /// * `coord` - The coordinate value to check against the intervals.
    ///
    /// # Returns
    ///
    /// A boolean value indicating whether the coordinate falls within any of the target intervals for the
    /// specified contig and strand.
    ///
    /// # Examples
    ///
    /// ```rust, ignore
    ///     ///
    /// let targets = Targets::new(TargetType::Direct(vec![
    ///     "Contig1,100,200,+".to_string(),
    ///     "Contig2,300,400,-".to_string(),
    /// ]));
    ///
    /// let record = CsvRecord {
    ///     contig: "Contig1".to_string(),
    ///     start: Some(150),
    ///     stop: Some(180),
    ///     strand: Some("+".to_string()),
    /// };
    ///
    /// let is_within_interval = record.get_coords("Contig1", "+", 160);
    ///
    /// assert!(is_within_interval);
    /// ```
    fn get_coords<T: ToString>(&self, contig: &str, strand: T, coord: usize) -> bool {
        let strand: Strand = strand.to_string().as_str().into();
        let intervals = self
            ._targets
            .get(&StrandWrapper(strand))
            .and_then(|inner_map| inner_map.get(contig));
        if let Some(intervals) = intervals {
            intervals
                .iter()
                .any(|&(start, end)| start <= coord && coord <= end)
        } else {
            false
        }
    }
}

impl Conf {
    /// Constructs a new [`Conf`] instance by parsing a TOML file.
    ///
    /// This function takes a TOML file path (`toml_path`) and reads its contents using [`std::fs::read_to_string`].
    /// The TOML content is then parsed into a `Table` using the `parse::<Table>` method. The [`Table`] represents
    /// the parsed TOML structure.
    ///
    /// The function initializes empty vectors `regions` and `barcodes` to hold the parsed regions and barcodes,
    /// respectively. It then checks if the parsed TOML structure contains the "regions" and "barcodes" sections.
    /// If the sections are present, the function iterates over the corresponding values and converts them into
    /// [`Region`] and [`Barcode`] structs, which are added to the `regions` and `barcodes` vectors, respectively.
    ///
    /// Finally, the function constructs and returns a new [`Conf`] instance with the populated `regions` and `barcodes`
    /// vectors. The `channels` field is set to 0, and the `_channel_map` field is initialized as an empty [`HashMap].
    ///
    /// # Arguments
    ///
    /// * `toml_path` - The path to the TOML file.
    ///
    /// # Returns
    ///
    /// A new [`Conf`] instance with the parsed regions and barcodes.
    ///
    /// # Panics
    ///
    /// This function panics if there is an error reading the TOML file or parsing its contents.
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// use my_module::Conf;
    ///
    /// let conf = Conf::new("config.toml");
    ///
    /// // Perform operations on the `conf` instance
    /// ```
    fn new(toml_path: impl AsRef<Path>) -> Conf {
        let toml_content = std::fs::read_to_string(toml_path).unwrap();
        let value = toml_content.parse::<Table>().unwrap();
        let mut regions = Vec::new();
        if let Some(parsed_regions) = value.get("regions") {
            let parsed_regions = parsed_regions.as_array().unwrap();
            for region in parsed_regions {
                let x = region.as_table().unwrap();
                let z: Region = Region {
                    condition: x.try_into().unwrap(),
                };
                regions.push(z);
            }
        }

        let mut barcodes = Vec::new();
        if let Some(parsed_barcodes) = value.get("barcodes") {
            let parsed_barcodes = parsed_barcodes.as_table().unwrap().iter();
            for (_, barcode_value) in parsed_barcodes {
                let x = barcode_value.as_table().unwrap();
                let z: Barcode = Barcode {
                    condition: x.try_into().unwrap(),
                };
                barcodes.push(z);
            }
        }
        Conf {
            channels: 0,
            regions,
            barcodes,
            _channel_map: HashMap::new(),
        }
    }
    // ToDo check for Unclassified and classified barcodes
}
/// Formats the sum of two numbers as string.
// #[pyfunction]
// fn sum_as_string(bam_path: String, toml_path: String) -> PyResult<String> {}

// /// A Python module implemented in Rust.
// #[pymodule]
// fn readfish_tools(_py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
//     Ok(())
// }

#[cfg(test)]
mod tests {
    // BEdfile, with not 6 rows, bedfile with wrong types, csv with wrong types, csv with more than 4 rws
    use toml::{Table, Value};

    use super::*;
    use std::fs;
    use std::path::PathBuf;

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
    fn test_get_csv_record_strand() {
        let record = CsvRecord {
            contig: "chr1".to_string(),
            start: Some(100),
            stop: Some(200),
            strand: Some("+".to_string()),
        };
        let strand = record.get_strand();
        assert_eq!(strand, Some(Strand::Forward));
        let record = CsvRecord {
            contig: "chr1".to_string(),
            start: Some(100),
            stop: Some(200),
            strand: Some("-1".to_string()),
        };
        let strand = record.get_strand();
        assert_eq!(strand, Some(Strand::Reverse))
    }

    #[test]
    fn test_insert_into_targets() {
        use std::collections::HashMap;
        let mut targets: HashMap<StrandWrapper, HashedTargets> = HashMap::new();
        let record = CsvRecord {
            contig: "chr1".to_string(),
            start: Some(100),
            stop: Some(200),
            strand: Some("+".to_string()),
        };
        Targets::insert_into_targets(&mut targets, &record, record.get_strand().unwrap());
        assert_eq!(targets.len(), 1);
        assert_eq!(
            targets.get(&StrandWrapper(Strand::Forward)).unwrap().len(),
            1
        );
        assert_eq!(
            targets
                .get(&StrandWrapper(Strand::Forward))
                .unwrap()
                .get("chr1")
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            targets
                .get(&StrandWrapper(Strand::Forward))
                .unwrap()
                .get("chr1")
                .unwrap()[0],
            (100, 200)
        );
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn read_toml() {
        let test_toml = get_test_file("RAPID_CNS2.toml");
        let toml_content = fs::read_to_string(test_toml).unwrap();
        let value = toml_content.parse::<Table>().unwrap();
        // println!("{:#?}", value);
        assert_eq!(
            value["regions"][0]["targets"].as_str(),
            Some("resources/panel_adaptive_nogenenames_20122021_hg38.bed")
        );
        assert!(match value["regions"][1]["targets"] {
            Value::Array(_) => true,
            Value::String(_) => false,
            _ => false,
        })
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn test_load_conf() {
        let test_toml = get_test_file("RAPID_CNS2.toml");
        let conf = Conf::new(test_toml);
        assert!(conf
            .regions
            .get(0)
            .map(|x| x.condition.name == "Rapid_CNS")
            .unwrap_or(false));
        assert!(conf
            .regions
            .get(1)
            .map(|x| x.condition.name == "Direct_CNS")
            .unwrap_or(false));
        assert!(conf
            .regions
            .get(1)
            .map(
                |x| x.condition.targets._targets[&StrandWrapper(Strand::Reverse)]["chr2"][0]
                    == (3000_usize, 4000_usize)
            )
            .unwrap_or(false));
        assert!(conf.barcodes.is_empty())
    }

    #[test]
    fn test_merge_intervals() {
        assert_eq!(
            Targets::_merge_intervals(&mut vec![
                (11, 15),
                (1, 3),
                (14, 17),
                (2, 4),
                (15, 100),
                (169, 173),
                (10, 29)
            ]),
            vec![(1, 4), (10, 100), (169, 173)]
        )
    }

    #[test]
    fn test_make_targets() {
        let targets: Targets = Targets::new(TargetType::Direct(vec![
            "chr1,10,20,+".to_string(),
            "chr1,15,30,+".to_string(),
        ]));
        assert_eq!(
            targets
                ._targets
                .get(&StrandWrapper(Strand::Forward))
                .unwrap()
                .get("chr1")
                .unwrap(),
            &vec![(10, 30)]
        )
    }

    #[test]

    fn test_get_coord() {
        let targets: Targets = Targets::new(TargetType::Direct(vec![
            "chr1,10,20,+".to_string(),
            "chr1,15,30,+".to_string(),
        ]));
        assert_eq!(
            targets
                ._targets
                .get(&StrandWrapper(Strand::Forward))
                .unwrap()
                .get("chr1")
                .unwrap(),
            &vec![(10, 30)]
        );
        assert!(targets.get_coords("chr1", Strand::Forward, 15));
        assert!(targets.get_coords("chr1", "+", 15));
        assert!(targets.get_coords("chr1", 1, 15));
        assert!(!targets.get_coords("chr1", 1, 40));
        assert!(!targets.get_coords("chr2", 1, 40));
        assert!(!targets.get_coords("chr1", "-", 15));
        assert!(!targets.get_coords("chr1", -1, 15));
    }

    #[test]
    fn test_get_coord_contig() {
        let targets: Targets = Targets::new(TargetType::Direct(vec!["chr1".to_string()]));
        assert_eq!(
            targets
                ._targets
                .get(&StrandWrapper(Strand::Forward))
                .unwrap()
                .get("chr1")
                .unwrap(),
            &vec![(0_usize, usize::MAX)]
        );
        assert!(targets.get_coords("chr1", Strand::Forward, 15));
        assert!(targets.get_coords("chr1", "+", 15));
        assert!(targets.get_coords("chr1", 1, 15));
        assert!(targets.get_coords("chr1", 1, 40));
        assert!(!targets.get_coords("chr2", 1, 40));
        assert!(targets.get_coords("chr1", "-", 15));
        assert!(targets.get_coords("chr1", -1, 15));
    }

    #[test]
    #[cfg_attr(miri, ignore)]
    fn load_barcoded_conf() {
        let test_toml = get_test_file("clockface.toml");
        let conf = Conf::new(test_toml);
        assert!(conf.regions.is_empty());
        assert!(conf
            .barcodes
            .get(0)
            .map(|x| x.condition.name == "barcode01")
            .unwrap_or(false));
        assert!(conf
            .barcodes
            .get(1)
            .map(|x| x.condition.name == "barcode02")
            .unwrap_or(false));
        assert!(conf
            .barcodes
            .get(2)
            .map(|x| x.condition.name == "barcode03")
            .unwrap_or(false));
        assert!(conf
            .barcodes
            .get(2)
            .map(
                |x| x.condition.targets._targets[&StrandWrapper(Strand::Reverse)]["NC_002516.2"][0]
                    == (0_usize, usize::MAX)
            )
            .unwrap_or(false))
    }
}
