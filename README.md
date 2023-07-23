# readfish-tools
A Python/rust wrapper for analysing the output of readfish runs. Can be used in conjunction with `readfish summarise`.
The overarching goal is to analyse a PAF, BAM or set of FASTQ files. This tool takes in a readfish TOML file, the sequencing summary output for the PAF file

Currently it is only possible to analyse PAF files, using the `readfish_tools.demultiplex_paf` function, whhich takes 3 parameters,
`toml_file`, `paf_file` and Optionally, `sequencing_summary`, which are file paths to the respective paths.
Currently if we do not find tags for the channel (ch) and optionally the barcode (ba) in the PAF tags, a sequencing summary file is required.

# Documentations
To build the rust documentation and view it -

```bash
cargo doc --no-deps --document-private-items --open
```
Python documentation is on the roadmap.

# Installing/Building
This should compile on X64 or Arm/arch64.

```bash
git clone https://github.com/Adoni5/readfish-tools
cd readfish-tools
# conda
mamba env install -f readfish_tools.yml
# pip
pip install -e .
```

# tests
To run rust integration, unit and doctests
```bash
cargo test
```

Python tests
```bash
pip install -e .[tests]
pytest --capture=tee-sys
```

# RoadMap
V0.0.2
    [#2](https://github.com/Adoni5/readfish-tools/issues/2)

More refined printing of summaries.
- [ ] More fields on the printout
- [ ] Comparison of given conditions
- [ ] Better grouping of stats
- [ ] Writing out of CSV files
- [ ] No stats options
- [ ] Take an iter of PAf records, rather than a full file.
- [ ] Python documentation


# Changelog
V0.0.1 - Basic printing out of stats, taken from custom tags in the PAF record, (ba, ch) or from a sequencing summary file.
