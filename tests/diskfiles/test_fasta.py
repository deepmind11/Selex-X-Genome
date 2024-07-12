from pathlib import Path

from diskfiles.countTables import CountTable


def test_build_count_table(my_Fasta_file):
    count_table = my_Fasta_file.build_count_table()
    assert isinstance(count_table, CountTable)
    assert count_table.file_path.exists()
    assert count_table.file_path == Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/tests/ENCODE/ENCFF476FQX.tsv.gz"
    )
