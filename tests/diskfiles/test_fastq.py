from pathlib import Path


def test_Fastq_subsample(my_Fastq_file):
    subsampled_fq = my_Fastq_file.subsample()
    assert subsampled_fq.file_path.exists()
    assert subsampled_fq.file_path == Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX_subsampled.fq"
    )
