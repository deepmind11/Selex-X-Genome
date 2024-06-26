from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from motifs.motif import Mononucleotide

import sqlite3
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from diskfiles.base import DiskFile


# Add fit porbound method for the count table
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
        script: Path = Path(__file__).parent.parent / Path("processCntTbl.sh"),
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

    # ! motif can be a general motif class
    def score(self, motif: Mononucleotide, search_tf: str):
        """Scores a motif against the count table and also generates the enrichment vs bin plot."""

        # Check if the table has already been scored. If so, just exit.
        with sqlite3.connect(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            # Get the count_table_id from the counttables table.
            try:
                cursor.execute(
                    """
                    SELECT "id" FROM "count_tables" WHERE "r1_file" IN (SELECT "id" FROM "files" WHERE "accession" = ?)
                    """,
                    (self.file_path.name[:11],),
                )
                count_table_id = cursor.fetchone()[0]
            except:
                raise ValueError("Could not find count table id in Database")

            try:
                cursor.execute(
                    """
                    SELECT COUNT(*) FROM "motif" WHERE "tf" = ? AND "search_tf" = ? AND "organism" = ? AND "count_table_id" = ?
                    """,
                    (
                        motif.tf,
                        search_tf,
                        motif.organism,
                        self.file_path.name[:11],
                    ),
                )
                count = cursor.fetchone()[0]
                if count != 0:
                    return
            except:
                raise ValueError("Failed to get count from Database")

        count_table_df = self.get_pandas_df()
        count_table_df["score"] = count_table_df["seq"].apply(motif.score_seq)
        count_table_df.sort_values(by="score", ascending=False, inplace=True)

        bin_df = CountTable.bin_count_table(count_table_df)

        # Plot Enrichment vs Bin
        # Step 2: Prepare the data
        x = list(range(bin_df.shape[0]))  # Example iterable for the x-axis
        y = list(bin_df["enrichment"])  # Example iterable for the y-axis

        # Step 3: Create the scatter plot
        fig, ax = plt.subplots()
        ax.scatter(x, y)

        # Step 4: Customize the plot
        ax.set_title(f"{motif.tf}_X_{search_tf} for {self.file_path.name[:11]}")
        ax.set_xlabel("Bin Number")
        ax.set_ylabel("Enrichment")

        fig.savefig(
            f"{str(self.file_path.parent)}/{self.file_path.name[:11]}_{motif.tf}_X_{search_tf}_Enr_vs_Bin.png"
        )

        # Plot Enrichment vs Score
        # Step 2: Prepare the data
        x = list(bin_df["avg_score"])  # Example iterable for the x-axis
        y = list(bin_df["enrichment"])  # Example iterable for the y-axis

        # Step 3: Create the scatter plot
        fig, ax = plt.subplots()
        ax.scatter(x, y)

        # Step 4: Customize the plot
        ax.set_title(f"{motif.tf}_X_{search_tf} for {self.file_path.name[:11]}")
        ax.set_xlabel("Score")
        ax.set_ylabel("Enrichment")

        fig.savefig(
            f"{str(self.file_path.parent)}/{self.file_path.name[:11]}_{motif.tf}_X_{search_tf}_Enr_vs_Score.png"
        )

        # Standard Scoring (PSAM not mean centered)
        count_table_df_st = self.get_pandas_df()
        count_table_df_st["score"] = count_table_df_st["seq"].apply(
            motif.score_seq_standard
        )
        count_table_df_st.sort_values(by="score", ascending=False, inplace=True)

        bin_df_st = CountTable.bin_count_table(count_table_df_st)

        # Plot Enrichment vs Bin
        # Step 2: Prepare the data
        x = list(range(bin_df_st.shape[0]))  # Example iterable for the x-axis
        y = list(bin_df_st["enrichment"])  # Example iterable for the y-axis

        # Step 3: Create the scatter plot
        fig, ax = plt.subplots()
        ax.scatter(x, y)

        # Step 4: Customize the plot
        ax.set_title(
            f"{motif.tf}_X_{search_tf} for {self.file_path.name[:11]} (Standard Scoring)"
        )
        ax.set_xlabel("Bin Number")
        ax.set_ylabel("Enrichment")

        fig.savefig(
            f"{str(self.file_path.parent)}/{self.file_path.name[:11]}_{motif.tf}_X_{search_tf}_Enr_vs_Bin_Standard.png"
        )

        # Plot Enrichment vs Score
        # Step 2: Prepare the data
        x = list(bin_df_st["avg_score"])  # Example iterable for the x-axis
        y = list(bin_df_st["enrichment"])  # Example iterable for the y-axis

        # Step 3: Create the scatter plot
        fig, ax = plt.subplots()
        ax.scatter(x, y)

        # Step 4: Customize the plot
        ax.set_title(f"{motif.tf}_X_{search_tf} for {self.file_path.name[:11]}")
        ax.set_xlabel("Score")
        ax.set_ylabel("Enrichment")

        fig.savefig(
            f"{str(self.file_path.parent)}/{self.file_path.name[:11]}_{motif.tf}_X_{search_tf}_Enr_vs_Score_Standard.png"
        )

        # Add the score to Database
        with sqlite3.connect(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            # Get the count_table_id from the counttables table.
            try:
                cursor.execute(
                    """
                    SELECT "id" FROM "count_tables" WHERE "r1_file" IN (SELECT "id" FROM "files" WHERE "accession" = ?)
                    """,
                    (self.file_path.name[:11],),
                )
                count_table_id = cursor.fetchone()[0]
            except:
                raise ValueError("Could not find file in Database")

            # Update DB
            type1 = "Mononucleotide"
            tf = motif.tf
            search_tf = search_tf
            organism = motif.organism
            count_table_id = count_table_id
            top_row = bin_df.iloc[0,].to_dict()
            top_row_st = bin_df_st.iloc[0,].to_dict()
            score = top_row["avg_score"]
            r0_count = top_row["r0_count"]
            r1_count = top_row["r1_count"]
            enrichment = top_row["enrichment"]
            score_st = top_row_st["avg_score"]
            r0_count_st = top_row_st["r0_count"]
            r1_count_st = top_row_st["r1_count"]
            enrichment_st = top_row_st["enrichment"]
            row = (
                type1,
                tf,
                search_tf,
                organism,
                count_table_id,
                score,
                r0_count,
                r1_count,
                enrichment,
                score_st,
                r0_count_st,
                r1_count_st,
                enrichment_st,
            )

            try:
                cursor.executemany(
                    """
                INSERT INTO "motif" ("type", "tf", "search_tf", "organism", "count_table_id", "score", "r0_count", "r1_count", "enrichment", "score_st", "r0_count_st", "r1_count_st", "enrichment_st") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                    [row],
                )
                conn.commit()
            except sqlite3.IntegrityError:
                pass
            except Exception as e:
                print("Error:", e)

        return bin_df.iloc[0,].to_dict()

    def get_probe_count(self):
        "Returns Total Probe Count, Probes in R0, and Probes in R1."
        df = self.get_pandas_df()
        return (df.shape[0], df[df["r0"] > 0].shape[0], df[df["r1"] > 0].shape[0])

    def update_database(self):
        "Adds information to the Selex_X_Genome.db database"
        # Connect to the database
        with sqlite3.connect(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            # Get the r1_file_id from the files table
            try:
                cursor.execute(
                    """
                    SELECT "id" FROM "files" WHERE "accession" = ?
                    """,
                    (self.file_path.name[:11],),
                )
                r1_file_id = cursor.fetchone()[0]
            except:
                raise ValueError("Could not find file in Database")

            probe_count, r0_count, r1_count = self.get_probe_count()
            # Value to be added to DB
            count_table = [(r1_file_id, probe_count, r0_count, r1_count)]
            # Update DB
            try:
                cursor.executemany(
                    """
                INSERT INTO "count_tables" ("r1_file", "probe_count", "r0_count", "r1_count") VALUES (?, ?, ?, ?)
                """,
                    count_table,
                )
                conn.commit()
            except sqlite3.IntegrityError:
                pass
            except Exception as e:
                print("Error:", e)

    def plot_enrichment_vs_bin(self, motif: Mononucleotide):
        count_table_df = self.get_pandas_df()
        count_table_df["score"] = count_table_df["seq"].apply(motif.score_seq)
        count_table_df.sort_values(by="score", ascending=False, inplace=True)

        bin_df = CountTable.bin_count_table(count_table_df)

        # Step 2: Prepare the data
        x = list(range(bin_df.shape[0]))  # Example iterable for the x-axis
        y = list(bin_df["enrichment"])  # Example iterable for the y-axis

        # Step 3: Create the scatter plot
        fig, ax = plt.subplots()
        ax.scatter(x, y)

        # Step 4: Customize the plot
        ax.set_title("Scatter Plot of Bin Number vs Enrichment")
        ax.set_xlabel("Bin Number")
        ax.set_ylabel("Enrichment")

        return fig, ax

    def fit_probound(self):
        pass
