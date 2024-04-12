import subprocess
from pathlib import Path

import pandas as pd
from diskfiles.base import DiskFile
from motifs.motif import Mononucleotide


class CountTable(DiskFile):
    """Class for CountTable (.tsv file) on disk."""

    def __init__(self, file_path: Path):
        super().__init__(file_path)

    @classmethod
    def create_from_fasta(cls, fasta1: Path, fasta2: Path, count_table_path: Path):
        if count_table_path.exists():
            return CountTable(count_table_path)
        elif (count_table_path.parent / Path(count_table_path.name + ".gz")).exists():
            return CountTable(
                count_table_path.parent / Path(count_table_path.name + ".gz")
            )

        """Creates a CountTable tsv file"""
        count_table_path.parent.mkdir(exist_ok=True, parents=True)
        count_table_path.touch(exist_ok=True)

        with open(count_table_path, "w") as output_file:
            subprocess.run(
                [
                    "/burg/home/hg2604/ProBound/pb",
                    "make-table",
                    str(fasta1),
                    str(fasta2),
                ],
                stdout=output_file,
            )

        return CountTable(count_table_path)

    def processTable(
        self,
        script: Path = Path(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/source/processCntTbl.sh"
        ),
    ):
        """Transforms count table in place, based on the script. Table should not be zipped(assert this). (In place)"""
        if not self.file_path.exists():
            raise FileNotFoundError(f"File {self.file_path} does not exist.")
        elif self.file_path.suffix == ".gz":
            return
        else:
            subprocess.run([str(script), str(self.file_path)])

    def get_pandas_df(self):
        """Returns a pandas dataframe of the count table"""
        if not self.file_path.exists():
            raise FileNotFoundError(f"File {self.file_path} does not exist.")
        else:
            return pd.read_csv(
                self.file_path,
                sep="\t",
                index_col=None,
                header=None,
                names=["seq", "r0", "r1"],
            )

    @staticmethod
    def bin_count_table(df, bin_size=1000):
        """df should be a count Table with score column"""  # Assert this
        i = 0
        avg_score = []
        r0_count = []
        r1_count = []
        while (
            i < df.shape[0] - 10
        ):  # Keeping a buffer, so that the last bin isn't super small
            bin = df.iloc[i : i + 1000, 1:].copy()
            avg_score.append(bin.iloc[:, 2].mean())
            r0_count.append(bin.iloc[:, 0].sum())
            r1_count.append(bin.iloc[:, 1].sum())
            i += 1000

        bin_df = pd.DataFrame(
            {"avg_score": avg_score, "r0_count": r0_count, "r1_count": r1_count}
        )
        bin_df["enrichment"] = bin_df["r1_count"] / bin_df["r0_count"]

        return bin_df

    def score(self, motif: Mononucleotide):
        """Scores a motif against the count table"""
        count_table_df = self.get_pandas_df()
        count_table_df["score"] = count_table_df["seq"].apply(motif.score_seq)
        count_table_df.sort_values(by="score", ascending=False, inplace=True)

        bin_df = CountTable.bin_count_table(count_table_df)

        return bin_df.iloc[0,].to_dict()
