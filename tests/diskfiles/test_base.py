from pathlib import Path

import pytest
from ENCODE.experiment import TFChipSeq


@pytest.fixture
def my_TFChipSeq():
    return TFChipSeq("ENCSR068HEE")


@pytest.fixture
def my_Exp_Library(my_TFChipSeq):
    my_Libraries = my_TFChipSeq.get_libraries()

    for library in my_Libraries:
        if library.accession == "ENCLB386VDZ":
            return library


@pytest.fixture
def my_File(my_Exp_Library):
    files = my_Exp_Library.get_Files()
    return files[0]


@pytest.fixture
def my_Fastq_file(my_File):
    fastq_file = my_File.download(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE"
    )
    return fastq_file


def test_Fastq_initialization(my_Fastq_file):
    assert not my_Fastq_file.deleted
    assert my_Fastq_file.zipped


def test_zipping(my_Fastq_file):

    my_Fastq_file.zip()  # Should do nothing?
    assert my_Fastq_file.zipped
    my_Fastq_file.unzip()
    assert not my_Fastq_file.zipped
    assert my_Fastq_file.file_path == Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX.fastq"
    )
    my_Fastq_file.zip()
    assert my_Fastq_file.zipped
    assert my_Fastq_file.file_path == Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX.fastq.gz"
    )
