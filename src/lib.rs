#![warn(missing_docs)]
#![allow(dead_code)]
//! # Readfish-tools
//!
//! `readfish-tools` is a collection of utilities to provide a standardised way of analysing
//! readfish runs that have been run. Currently the accepted analysable inputs are sequencing summary files,
//! BAM of all produced FASTQ, and the `TOML` file that was used in the readfish run.
//!
//! #

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

#[derive(Debug)]
enum Action {
    Unblock,
    StopReceiving,
    Proceed,
}
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

/// The Region struct holds the settings lifted from the TOML file, for each
/// region of the flowcell.
#[derive(Debug)]
struct _Condition {
    /// The name of the region
    name: String,
    min_chunks: u8,
    max_chunks: u8,
    targets: Targets,
    single_off: Action,
    single_on: Action,
    multi_off: Action,
    multi_on: Action,
    no_map: Action,
    no_seq: Action,
}

#[derive(Debug)]
struct Region {
    condition: _Condition,
}

#[derive(Debug)]
struct Barcode {
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

#[derive(Debug, Hash, PartialEq)]
enum Strand {
    Forward,
    Reverse,
}

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
#[derive(Clone, Debug)]
enum TargetType {
    Direct(Vec<String>),
    ViaFile(PathBuf),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct BedRecord {
    contig: String,
    start: usize,
    stop: usize,
    _name: String,
    _score: String,
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
    /// Returns true if there were start and stop coordinates including in this CSV records.
    /// Example
    /// ```rust
    /// # #[macro_use] extern crate readfish_tools;
    /// let csv: readfish_tools::CsvRecord  = readfish_tools::CsvRecord{contig: "contig".to_string(), start: Some(10), stop: Some(20), strand: Some("+".to_string())};
    /// assert!(csv.has_coords())
    /// ```
    pub fn has_coords(&self) -> bool {
        self.start.is_some() & self.stop.is_some()
    }

    fn get_coords(&self) -> (usize, usize) {
        if self.has_coords() {
            (self.start.unwrap(), self.stop.unwrap())
        } else {
            (0, usize::MAX)
        }
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

#[derive(Debug)]
struct Conf {
    channels: usize,
    regions: Vec<Region>,
    barcodes: Vec<Barcode>,
    _channel_map: HashMap<u16, u8>,
}

#[derive(Debug)]
struct Targets {
    value: TargetType,
    _targets: HashMap<StrandWrapper, HashMap<String, Vec<(usize, usize)>>>,
}

impl Targets {
    fn new(targets: TargetType) -> Targets {
        let t = targets.clone();
        Targets {
            value: targets,
            _targets: Targets::from_parsed_toml(t),
        }
    }

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
                            record.strand.as_ref().unwrap().as_str().into(),
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
                            record.strand.as_ref().unwrap().as_str().into(),
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

    /// For given genomic coordinate return true if on target, false if off target.
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
    fn load_conf() {
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

    #[test]
    fn get_channels() {
        println!("{:#?}", channels::MINION_CHANNELS[&1]);
    }
}
