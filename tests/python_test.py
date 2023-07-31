from pathlib import Path
import copy
import re

import pytest
from mappy import fastx_read
import mappy_rs
from readfish_tools import summarise_paf, ReadfishSummary


RESOURCES = Path(__file__).parent.resolve().parent.resolve() / "resources/"
TOML_FILE = RESOURCES / "human_barcode.toml"
PAF_FILE = RESOURCES / "test_paf_barcode05_NA12878.chr.paf"
SEQ_SUM_FILE = RESOURCES / "seq_sum_PAK09329.txt"
MMI_FILE = "/home/adoni5/Documents/Bioinformatics/refs/hg38_no_alts_22.mmi"


def get_fq(directory):
    """
    Given a directory, return a generator of fastq files.

    Parameters:
        directory (str or Path): The directory path to search for fastq files.

    Yields:
        str: A path to a fastq file found in the given directory or its subdirectories.

    Examples:
        >>> for file_path in get_fq("resouces"):
        ...     print(file_path)
        /path/to/directory/sample1.fastq
        /path/to/directory/sample2.fastq.gz
        ...

    Note:
        The function searches for files with extensions .fastq, .fastq.gz, .fq, .fq.gz
        in the specified directory and its subdirectories.
    """
    types = {".fastq", ".gz", ".fq"}
    files = (
        str(p.resolve())
        for p in Path(directory).glob("**/*")
        if set(map(str.lower, p.suffixes)).intersection(types)
    )
    yield from files


@pytest.fixture
def toml_file_path():
    return TOML_FILE


@pytest.fixture
def paf_file_path():
    return PAF_FILE


@pytest.fixture
def seq_sum_file_path():
    return SEQ_SUM_FILE


@pytest.fixture
def toml_file_str():
    return str(TOML_FILE)


@pytest.fixture
def paf_file_str():
    return str(PAF_FILE)


@pytest.fixture
def seq_sum_file_str():
    return str(SEQ_SUM_FILE)


@pytest.fixture
def mmi_file():
    return str(MMI_FILE)


@pytest.fixture
def al(mmi_file):
    return mappy_rs.Aligner(mmi_file)


def _prep_fastq():
    """
    Function extracts key-value pairs
    from the comments present in the fastq files and yields a tuple for each file with the processed sequence,
    name, channel, and barcode information.

    Yields
    ------
    tuple
        A tuple containing a dictionary with the processed sequence and a 3-tuple with the following information:
        1. The name extracted from the fastq file.
        2. The channel (integer) extracted from the comments in the fastq file.
        3. The barcode (optional, None if not present) extracted from the comments in the fastq file.

    Example
    -------
    .. code-block:: python

        for seq_info, metadata in _prep_fastq():
            processed_sequence = seq_info["seq"]
            name, channel, barcode = metadata
            print(f"Name: {name}, Channel: {channel}, Barcode: {barcode}, Sequence: {processed_sequence}")

    Output:
    .. code-block:: python

        Name: some_name, Channel: 123, Barcode: barcode1, Sequence: ACGT
        Name: another_name, Channel: 456, Barcode: None, Sequence: ACGT
    """
    # Define a regex pattern to capture each side of the "=" sign

    pattern = r"(\w+)=([^ =]+)"
    pattern = re.compile(pattern)

    for file in get_fq("resources/barcoded_fastq"):
        for name, seq, _qual, comment in fastx_read(file, read_comment=True):
            # Find all matches of the pattern in the input string
            comments = dict(pattern.findall(comment))
            channel = int(comments["ch"])
            barcode = comments.get("barcode", None)
            yield {
                "seq": seq,
                "name": name,
                "channel": channel,
                "barcode": barcode,
                "seq_len": len(seq),
            }


def tupleise(input: dict):
    return (input["name"], input["channel"], input.get("barcode", None))


def yield_alignments(al):
    for mappings, _input in al.map_batch(_prep_fastq()):
        for mapping in mappings:
            metadata = tupleise(_input)
            seq_len = _input.get("seq_len")
            yield (f"{metadata[0]}\t{seq_len}\t{mapping})", tupleise(_input))


def test_summarise_class():
    rfs = ReadfishSummary()


def test_map_and_parse(al):
    """
    Test mapping and parsing using the summarise class API of readfish tools
    """
    al.enable_threading(2)

    rfs = ReadfishSummary()
    rfs.with_toml_conf(TOML_FILE)
    rfs.parse_paf_from_iter(yield_alignments(al))
    rfs.print_summary()


def test_get_fq(tmpdir):
    # Create temporary directory for testing
    directory = tmpdir.mkdir("test_dir")

    # Create some sample fastq files with different extensions
    directory.join("sample1.fastq").write("Sample content")
    directory.join("sample2.fastq.gz").write("Sample content")
    directory.join("sample3.fq").write("Sample content")
    directory.join("sample4.fq.gz").write("Sample content")
    directory.join("sample5.txt").write("Not a fastq file")

    # Get the list of fastq files using the function being tested
    fastq_files = list(get_fq(directory))

    # Assert that the function returns the correct number of files
    assert len(fastq_files) == 4

    # Assert that all the expected fastq files are in the list
    assert str(directory.join("sample1.fastq")) in fastq_files
    assert str(directory.join("sample2.fastq.gz")) in fastq_files
    assert str(directory.join("sample3.fq")) in fastq_files
    assert str(directory.join("sample4.fq.gz")) in fastq_files

    # Assert that a non-fastq file is not included in the list
    assert str(directory.join("sample5.txt")) not in fastq_files


def test_demultiplex_pathlib(toml_file_path, paf_file_path, seq_sum_file_path):
    summarise_paf(
        toml_file_path,
        paf_file_path,
        seq_sum_file_path,
    )
