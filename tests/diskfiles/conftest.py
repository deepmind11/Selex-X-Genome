from pathlib import Path

import pytest
from motifs.parse_motifcentral_json import MOTIFCENTRAL


@pytest.fixture
def my_TFChipSeq():
    from ENCODE.experiment import TFChipSeq

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


@pytest.fixture
def my_subsampled_Fastq_file(my_Fastq_file):
    subsampled_fq = my_Fastq_file.subsample()
    return subsampled_fq


@pytest.fixture
def my_Fasta_file(my_subsampled_Fastq_file):
    fasta = my_subsampled_Fastq_file.transform_to_fasta()
    return fasta


@pytest.fixture
def my_countTable():
    from diskfiles.countTables import CountTable

    return CountTable(
        Path(
            "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX.tsv.gz"
        )
    )


@pytest.fixture
def my_countTable_sample():
    from diskfiles.countTables import CountTable

    return CountTable(
        Path(
            "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX_sample.tsv.gz"
        )
    )


@pytest.fixture
def my_mono_motif():
    from motifs.motif import Mononucleotide

    return Mononucleotide.create_from_motif_central(MOTIFCENTRAL[222])
