# Experiment:ENCSR068HEE  Control:ENCSR608IVH
import sqlite3
from pathlib import Path

import pytest
from ENCODE.base import RunType
from ENCODE.experiment import TFChipSeq
from ENCODE.files import SE_File


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
def my_Control_Library(my_TFChipSeq):
    my_Control = my_TFChipSeq.get_controls()
    control_libraries = my_Control.get_libraries()

    for library in control_libraries:
        if library.accession == "ENCLB312FJK":
            return library


def test_Exp_Library_Initialization(my_Exp_Library):
    assert my_Exp_Library.accession == "ENCLB386VDZ"
    assert my_Exp_Library.antibody == "ENCAB719MQZ"
    assert my_Exp_Library.technical_replicate_number == 2
    assert my_Exp_Library.biological_replicate_number == 1
    assert my_Exp_Library.experiment == "ENCSR068HEE"
    assert my_Exp_Library.biosample == "ENCBS672IKS"
    assert my_Exp_Library.run_type == RunType.SE
    assert my_Exp_Library.files is not None
    assert not my_Exp_Library.control


def test_Control_Library_Initialization(my_Control_Library):
    assert my_Control_Library.accession == "ENCLB312FJK"
    assert my_Control_Library.antibody is None
    assert my_Control_Library.technical_replicate_number == 2
    assert my_Control_Library.biological_replicate_number == 2
    assert my_Control_Library.experiment == "ENCSR608IVH"
    assert my_Control_Library.biosample == "ENCBS672IKS"
    assert my_Control_Library.run_type == RunType.SE
    assert my_Control_Library.files is not None
    assert my_Control_Library.control


def test_Exp_Library_get_files(my_Exp_Library):
    files = my_Exp_Library.get_Files()  # All fastq
    assert isinstance(files, list)
    assert len(files) == 1
    for file in files:
        assert isinstance(file, SE_File)
        assert file.file_format == "fastq"
        assert file.library == "ENCLB386VDZ"
        assert file.run_type == RunType.SE
        assert file.paired_end is None
        assert file.paired_with is None
        assert not file.control


@pytest.fixture
def my_File(my_Exp_Library):
    files = my_Exp_Library.get_Files()
    return files[0]


@pytest.fixture
def my_Control_File(my_Control_Library):
    files = my_Control_Library.get_Files()
    for file in files:
        if file.accession == "ENCFF392SII":
            return file


def test_File_Initialization(my_File):
    assert my_File.file_format == "fastq"
    assert my_File.library == "ENCLB386VDZ"
    assert my_File.run_type == RunType.SE
    assert my_File.paired_end is None
    assert my_File.paired_with is None
    assert my_File.accession == "ENCFF476FQX"
    assert my_File.read_length == 76
    assert my_File.experiment == "ENCSR068HEE"
    assert not my_File.control
    assert my_File.biosample == "ENCBS672IKS"  # Sample from which reads were derived
    assert my_File.technical_replicate_number == 2
    assert my_File.biological_replicate_number == 1
    assert my_File.antibody == "ENCAB719MQZ"


def test_File_download(my_File):
    fastq_file = my_File.download(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE"
    )
    assert fastq_file.file_path.exists()
    assert fastq_file.file_path == Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX.fastq.gz"
    )

    assert fastq_file.accession == "ENCFF476FQX"
    assert fastq_file.platform == "Illumina NextSeq 500"
    assert fastq_file.read_length == 76
    assert fastq_file.experiment == "ENCSR068HEE"
    assert fastq_file.library == "ENCLB386VDZ"
    assert fastq_file.biosample == "ENCBS672IKS"
    assert fastq_file.technical_replicate_number == 2
    assert fastq_file.biological_replicate_number == 1
    assert not fastq_file.control
    assert fastq_file.antibody == "ENCAB719MQZ"
    assert fastq_file.run_type == RunType.SE  # This is either SE or PE
    assert fastq_file.paired_end is None  # Will be none for SE (Assert this)
    assert fastq_file.paired_with is None


def test_Control_File_download(my_Control_File):
    fastq_file = my_Control_File.download(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE"
    )
    assert fastq_file.file_path.exists()
    assert fastq_file.file_path == Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF392SII.fastq.gz"
    )

    assert fastq_file.accession == "ENCFF392SII"
    assert fastq_file.platform == "Illumina NextSeq 500"
    assert fastq_file.read_length == 76
    assert fastq_file.experiment == "ENCSR608IVH"
    assert fastq_file.library == "ENCLB312FJK"
    assert fastq_file.biosample == "ENCBS672IKS"
    assert fastq_file.technical_replicate_number == 2
    assert fastq_file.biological_replicate_number == 2
    assert fastq_file.control
    assert fastq_file.antibody is None
    assert fastq_file.run_type == RunType.SE  # This is either SE or PE
    assert fastq_file.paired_end is None  # Will be none for SE (Assert this)
    assert fastq_file.paired_with is None


def test_update_db_library(my_Exp_Library):
    my_Exp_Library.update_database()

    with sqlite3.connect(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
    ) as conn:

        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT "antibody", "biosample", "technical_rep_number", "biological_rep_number" FROM "libraries" WHERE "accession" = ?
            """,
            (my_Exp_Library.accession,),
        )
        results = cursor.fetchone()
        antibody, biosample, technical_rep_number, biological_rep_number = results

    assert antibody == my_Exp_Library.antibody
    assert biosample == my_Exp_Library.biosample
    assert technical_rep_number == my_Exp_Library.technical_replicate_number
    assert biological_rep_number == my_Exp_Library.biological_replicate_number


def test_update_db_file(my_File):
    my_File.update_database()

    with sqlite3.connect(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
    ) as conn:

        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT "read_length", "paired_end" FROM "files" WHERE "accession" =?
            """,
            (my_File.accession,),
        )
        results = cursor.fetchone()
        read_length, paired_end = results

    assert read_length == my_File.read_length
    assert paired_end == "n"
