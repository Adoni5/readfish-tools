from pathlib import Path
import copy

import pytest

RESOURCES = Path(__file__).parent.resolve().parent.resolve() / "resources/"
TOML_FILE = RESOURCES / "human_barcode.toml"
PAF_FILE = RESOURCES / "test_paf_barcode05_NA12878.chr.paf"
SEQ_SUM_FILE = RESOURCES / "seq_sum_PAK09329.txt"

from readfish_tools import demultiplex_paf


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


def test_demultiplex_pathlib(toml_file_path, paf_file_path, seq_sum_file_path):
    demultiplex_paf(
        toml_file_path,
        paf_file_path,
        seq_sum_file_path,
    )
