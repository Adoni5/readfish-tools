# readfish-tools
A Python/rust wrapper for analysing the output of readfish runs. Can be used in conjunction with `readfish summarise`.
The overarching goal is to analyse a PAF, BAM or set of FASTQ files. This tool takes in a readfish TOML file, the sequencing summary output for the PAF file

Currently it is only possible to analyse PAF files, using the `readfish_tools.demultiplex_paf` function or the `ReadfishSummary` class. Usage for both is described below.

# Documentations
To build the rust documentation and view it -

```bash
cargo doc --no-deps --document-private-items --open
```
Python documentation is on the roadmap.

# Installing/Building
This should compile on X64 or Arm/arch64.

:::warning
⚡ Note that if installed with test dependencies, (which is the default in the conda env yaml) `mappy_rs` will be installed, which is NOT Arm/aarch64 compatible (yet).
:::

```bash
git clone https://github.com/Adoni5/readfish-tools
cd readfish-tools
# conda (HAS test dependencies)
mamba env install -f readfish_tools_test.yml
# or via pip, without test dependencies
pip install -e .
```

# Usage
Can be imported either as a Summary class which can be worked with, or a one shot function, which consumes all records in a given file. There are different limitations to each approach.

## Summary class
```python
from readfish_tools import ReadfishSummary

rfs = ReadfishSummary()
rfs.with_toml_conf(<TOML_FILE>)
rfs.parse_paf_from_iter(<iterable of tuple of (pafline, (read_id, channel number, Optional[barcode name]))>)
rfs.print_summary()
```

First it is necessary to initialise the class. The class has methods to set the configuration TOML file. There is also a method to set a path to a sequencing summary file, which is currently **unimplemented**.

As such the only way to successfully use the class of this moment is by ensuring that the iterator provided provides the additional tuple that contains the metadata, as the second element of the tuple, with the paf record as the first element.

For example a valid tuple iterator could look like:

```python
iter([("read123  100 0   100 +   contig123   300 0   300 200 200 50  ch=1", ("read123", 1, None))])
# or if barcoded
iter([("read123  100 0   100 +   contig123   300 0   300 200 200 50  ch=1", ("read123", 1, "barcode01"))])

```

It is possible to call the `parse_paf_from_iter` method multiple times, to parse multiple files, or to parse a single file in chunks. It is also possible to call the `print_summary` method more than once, and
the summary printed will represent the given parsed data at any point when called.

`print_summary` prints to stdout, and will print a table created by the `prettytable.rs` crate.

On the roadmap is a function to return manipulatable `ConditionSummary` and `ContigSumary` classes which can be manipulated in python.

## One shot function
```python
from readfish_tools import summarise_paf
summarise_paf(<TOML_PATH>, <PAF_FILE_PATH>, <SEQUENCING_SUMMARY_PATH>)
# Summarised table
#+---------------------------+-------------+----------------+--------------+-------------+------------+-----------+-----------+-----------+------------+
#| Condition                 | Total reads | # Off-target   | # On-target  | Total Yield | Off Target | On Target | Mean read | On target | Off target |
#|                           |             | reads          | reads        |             |  Yield     |  yield    |  length   | Mean read | Mean read  |
#|                           |             |                |              |             |            |           |           |  length   |  length    |
#+---------------------------+-------------+----------------+--------------+-------------+------------+-----------+-----------+-----------+------------+
#| barcode05_NA12878_tst-170 | 4,236       | 4,210 (99.39%) | 26 (0.61%)   | 3.90 Mb     | 3.79 Mb    | 111.62 Kb | 969 b     | 4.29 Kb   | 885 b      |
#+---------------------------+-------------+----------------+--------------+-------------+------------+-----------+-----------+-----------+------------+
#+----------------+---------------+-------------+-----------+-------------+-----------+--------------+-----------+-----------+------------+
#| Condition Name | barcode05_NA12878_tst-170   |           |             |           |              |           |           |            |
#+----------------+---------------+-------------+-----------+-------------+-----------+--------------+-----------+-----------+------------+
#| Contig         | Contig Length | Read count  | Yield     | Mean        | On Target | Off          | Mean read | On target | Off target |
#|                |               |             |           | Read Length |  Reads    | Target Reads |  length   | Mean read | Mean read  |
#|                |               |             |           |             |           |              |           |  length   |  length    |
#+----------------+---------------+-------------+-----------+-------------+-----------+--------------+-----------+-----------+------------+
#| chr1           | 248,956,422   | 352         | 335.21 Kb | 944 b       | 0         | 352          | 944 b     | 0 b       | 944 b      |
#+----------------+---------------+-------------+-----------+-------------+-----------+--------------+-----------+-----------+------------+
#...
```
The `summarise_paf` function takes 3 parameters, `toml_file`, `paf_file` and Optionally, `sequencing_summary`, which are file paths to the respective paths.
Currently if we do not find **custom** tags for the channel (ch) and optionally the barcode (ba) in the PAF tags, a sequencing summary file is required.
### Limitations

Currently, if a sequencing summary file is provided, a record buffer of 100,000 rows is filled. If the Paf record being analysed is not found in this buffer, the buffer rolls along the file, removing the oldest line when a new line is read.
Therefore, if the PAF file being analysed is not in the order in which reads were base-called (with 100,000 reads leeway), the analysis will not work properly, with some reads being skipped.
This is most likely to be a problem on barcoded runs.

# tests
To run rust integration, unit and doctests
```bash
cargo test
```

Python tests
```bash
pip install -e .[tests]
pytest -sv
```

# RoadMap
V0.0.2
    [#2](https://github.com/Adoni5/readfish-tools/issues/2)

More refined printing of summaries.
- [ ] More fields on the printout
- [ ] Comparison of given conditions
- [X] Better grouping of stats
- [ ] Writing out of CSV files
- [ ] No stats options
- [X] Take an iter of PAf records, rather than a full file.
- [ ] Python documentation


# Changelog
V0.0.1 - Basic printing out of stats, taken from custom tags in the PAF record, (ba, ch) or from a sequencing summary file.
