use readfish_tools::demultiplex_paf;

// importing the common code for tests.
mod common;

#[test]
fn test_region_based_paf_demultiplex() {
    // using common code.
    let paf = common::get_test_file("test_paf_barcode05_NA12878.paf");
    let seq_sum = common::get_test_file("seq_sum_PAK09329.txt")
        .as_os_str()
        .to_str()
        .unwrap()
        .to_string();
    let toml_path = common::get_test_file("RAPID_CNS2.toml");
    demultiplex_paf(toml_path, paf, Some(seq_sum))
}
