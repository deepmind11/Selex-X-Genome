from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from motifs.motif import Mononucleotide

import sqlite3
import subprocess
from collections import defaultdict
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
        """
        Scores a motif against the count table using the fit_id and ProboundTools.
        Also generates the enrichment vs bin plot and enrichment vs score plots.
        """

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

        # unzip count_table
        self.unzip()
        # score it
        subprocess.run(
            [
                str(Path(__file__).parent.parent / Path("scoreCntTbl.sh")),
                str(self.file_path),
                motif.fit_id,
            ]
        )
        # zip it again
        self.zip()

        count_table_df = self.get_pandas_df()

        score_series = pd.read_csv(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_{motif.fit_id}.txt"),
            header=None,
        )

        score_series = score_series.iloc[:, 0]

        count_table_df["score"] = score_series
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
            score = top_row["avg_score"]
            r0_count = top_row["r0_count"]
            r1_count = top_row["r1_count"]
            enrichment = top_row["enrichment"]
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
            )

            try:
                cursor.executemany(
                    """
                INSERT INTO "motif" ("type", "tf", "search_tf", "organism", "count_table_id", "score", "r0_count", "r1_count", "enrichment") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
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

    def plot_enrichment_vs_bin(self, motif_id: str):
        count_table_df = self.get_pandas_df()

        score_series = pd.read_csv(
            self.file_path.parent / Path(f"{self.file_path.name[:11]}_{motif_id}.txt"),
            header=None,
        )

        score_series = score_series.iloc[:, 0]

        count_table_df["score"] = score_series
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

        fig.savefig(
            f"{str(self.file_path.parent)}/{self.file_path.name[:11]}_{motif_id}_X_{motif_id}_Enr_vs_Bin.png"
        )

        return fig, ax

    def plot_enrichment_vs_score(self, motif_id: str):

        count_table_df = self.get_pandas_df()
        score_series = pd.read_csv(
            self.file_path.parent / Path(f"{self.file_path.name[:11]}_{motif_id}.txt"),
            header=None,
        )

        score_series = score_series.iloc[:, 0]

        count_table_df["score"] = score_series
        count_table_df.sort_values(by="score", ascending=False, inplace=True)

        bin_df = CountTable.bin_count_table(count_table_df)

        # Plot Enrichment vs Score
        # Step 2: Prepare the data
        x = list(bin_df["avg_score"])  # Example iterable for the x-axis
        y = list(bin_df["enrichment"])  # Example iterable for the y-axis

        # Step 3: Create the scatter plot
        fig, ax = plt.subplots()
        ax.scatter(x, y)

        # Step 4: Customize the plot
        ax.set_title(f"{motif_id}_X_{motif_id} for {self.file_path.name[:11]}")
        ax.set_xlabel("Score")
        ax.set_ylabel("Enrichment")

        fig.savefig(
            f"{str(self.file_path.parent)}/{self.file_path.name[:11]}_{motif_id}_X_{motif_id}_Enr_vs_Score.png"
        )

        return fig, ax

    def get_mn_base_composition(self):
        "Saves the base composition of each sequence as a tsv file. [A,C,G,T]"
        table_df = self.get_pandas_df()

        base_counts = {"A": [], "C": [], "G": [], "T": []}

        for index, row in table_df.iterrows():
            seq = row["seq"]
            probe_length = len(seq)
            count = defaultdict(int)
            for base in seq:
                count[base] += 1

            base_counts["A"].append(count["A"] / probe_length)
            base_counts["C"].append(count["C"] / probe_length)
            base_counts["G"].append(count["G"] / probe_length)
            base_counts["T"].append(count["T"] / probe_length)

        base_comp = pd.DataFrame.from_dict(base_counts)
        base_comp.to_csv(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_mn_composition.tsv"),
            sep="\t",
        )

        return base_comp

    def get_dn_base_composition(self):
        "Saves the dinucleotide base composition of each sequence as a tsv file. [A,C,G,T]"

        table_df = self.get_pandas_df()

        base_counts = {
            "AA": [],
            "AC": [],
            "AG": [],
            "AT": [],
            "CA": [],
            "CC": [],
            "CG": [],
            "CT": [],
            "GA": [],
            "GC": [],
            "GG": [],
            "GT": [],
            "TA": [],
            "TC": [],
            "TG": [],
            "TT": [],
        }

        for index, row in table_df.iterrows():
            seq = row["seq"]
            probe_length = len(seq)
            count = defaultdict(int)
            for index, _ in enumerate(seq):
                if index == len(seq) - 1:
                    break
                count[seq[index : index + 2]] += 1

            base_counts["AA"].append(count["AA"] / probe_length)
            base_counts["AC"].append(count["AC"] / probe_length)
            base_counts["AG"].append(count["AG"] / probe_length)
            base_counts["AT"].append(count["AT"] / probe_length)
            base_counts["CA"].append(count["CA"] / probe_length)
            base_counts["CC"].append(count["CC"] / probe_length)
            base_counts["CG"].append(count["CG"] / probe_length)
            base_counts["CT"].append(count["CT"] / probe_length)
            base_counts["GA"].append(count["GA"] / probe_length)
            base_counts["GC"].append(count["GC"] / probe_length)
            base_counts["GG"].append(count["GG"] / probe_length)
            base_counts["GT"].append(count["GT"] / probe_length)
            base_counts["TA"].append(count["TA"] / probe_length)
            base_counts["TC"].append(count["TC"] / probe_length)
            base_counts["TG"].append(count["TG"] / probe_length)
            base_counts["TT"].append(count["TT"] / probe_length)

        base_comp = pd.DataFrame.from_dict(base_counts)
        base_comp.to_csv(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_dn_composition.tsv"),
            sep="\t",
        )

        return base_comp

    def plot_perbin_mn_composition(self, motif_id: str):

        table_df = self.get_pandas_df()

        score_series = pd.read_csv(
            self.file_path.parent / Path(f"{self.file_path.name[:11]}_{motif_id}.txt"),
            header=None,
        )

        score_series = score_series.iloc[:, 0]

        table_df["score"] = score_series

        mn_comp = pd.read_csv(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_mn_composition.tsv"),
            sep="\t",
            index_col=False,
        )

        # Concatenating the Columns

        merged_df = pd.concat([table_df, mn_comp], axis=1, ignore_index=False)
        merged_df.sort_values(by="score", ascending=False, inplace=True)
        merged_df.reset_index(drop=True, inplace=True)

        # Start Binning

        i = 0
        avg_score = []
        r0_count = []
        r1_count = []
        A = []
        C = []
        G = []
        T = []

        while (
            i < merged_df.shape[0] - 10
        ):  # Keeping a buffer, so that the last bin isn't super small
            bin = merged_df.iloc[i : i + 1000,].copy()

            avg_score.append(bin.iloc[:, 2].mean())
            r0_count.append(bin.iloc[:, 0].sum())
            r1_count.append(bin.iloc[:, 1].sum())
            A.append(bin.iloc[:, 5].mean())
            C.append(bin.iloc[:, 6].mean())
            G.append(bin.iloc[:, 7].mean())
            T.append(bin.iloc[:, 8].mean())
            i += 1000

        # Make the plot
        x_axis = range(len(A))

        # plot
        fig, ax = plt.subplots()

        ax.set_title(
            f"Mono-Nucleotide composition across bins for {self.file_path.name[:11]}"
        )

        ax.plot(x_axis, A, color="b", label="A")
        ax.plot(x_axis, C, color="g", label="C")
        ax.plot(x_axis, G, color="r", label="G")
        ax.plot(x_axis, T, color="y", label="T")
        ax.legend()

        plt.savefig(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_{motif_id}_mn_composition.png")
        )

        return

    def plot_perbin_dn_composition(self, motif_id):

        table_df = self.get_pandas_df()

        score_series = pd.read_csv(
            self.file_path.parent / Path(f"{self.file_path.name[:11]}_{motif_id}.txt"),
            header=None,
        )

        score_series = score_series.iloc[:, 0]

        table_df["score"] = score_series

        dn_comp = pd.read_csv(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_dn_composition.tsv"),
            sep="\t",
            index_col=False,
        )

        # Concatenating the Columns

        merged_df = pd.concat([table_df, dn_comp], axis=1, ignore_index=False)
        merged_df.sort_values(by="score", ascending=False, inplace=True)
        merged_df.reset_index(drop=True, inplace=True)

        # Start Binning

        i = 0
        avg_score = []
        r0_count = []
        r1_count = []
        AA = []
        AC = []
        AG = []
        AT = []
        CA = []
        CC = []
        CG = []
        CT = []
        GA = []
        GC = []
        GG = []
        GT = []
        TA = []
        TC = []
        TG = []
        TT = []

        while (
            i < merged_df.shape[0] - 10
        ):  # Keeping a buffer, so that the last bin isn't super small
            bin = merged_df.iloc[i : i + 1000,].copy()

            avg_score.append(bin.iloc[:, 2].mean())
            r0_count.append(bin.iloc[:, 0].sum())
            r1_count.append(bin.iloc[:, 1].sum())
            AA.append(bin.iloc[:, 5].mean())
            AC.append(bin.iloc[:, 6].mean())
            AG.append(bin.iloc[:, 7].mean())
            AT.append(bin.iloc[:, 8].mean())
            CA.append(bin.iloc[:, 9].mean())
            CC.append(bin.iloc[:, 10].mean())
            CG.append(bin.iloc[:, 11].mean())
            CT.append(bin.iloc[:, 12].mean())
            GA.append(bin.iloc[:, 13].mean())
            GC.append(bin.iloc[:, 14].mean())
            GG.append(bin.iloc[:, 15].mean())
            GT.append(bin.iloc[:, 16].mean())
            TA.append(bin.iloc[:, 17].mean())
            TC.append(bin.iloc[:, 18].mean())
            TG.append(bin.iloc[:, 19].mean())
            TT.append(bin.iloc[:, 20].mean())
            i += 1000

        # Make the plot
        x_axis = range(len(AA))

        # plot
        fig, ax = plt.subplots()

        # Color map
        cmap = plt.get_cmap("tab20")

        ax.set_title(
            f"Di-Nucleotide composition across bins for {self.file_path.name[:11]}"
        )

        ax.plot(x_axis, AA, color=cmap(0), label="AA")
        ax.plot(x_axis, AC, color=cmap(1), label="AC")
        ax.plot(x_axis, AG, color=cmap(2), label="AG")
        ax.plot(x_axis, AT, color=cmap(3), label="AT")
        ax.plot(x_axis, CA, color=cmap(4), label="CA")
        ax.plot(x_axis, CC, color=cmap(5), label="CC")
        ax.plot(x_axis, CG, color=cmap(6), label="CG")
        ax.plot(x_axis, CT, color=cmap(7), label="CT")
        ax.plot(x_axis, GA, color=cmap(8), label="GA")
        ax.plot(x_axis, GC, color=cmap(9), label="GC")
        ax.plot(x_axis, GG, color=cmap(10), label="GG")
        ax.plot(x_axis, GT, color=cmap(11), label="GT")
        ax.plot(x_axis, TA, color=cmap(12), label="TA")
        ax.plot(x_axis, TC, color=cmap(13), label="TC")
        ax.plot(x_axis, TG, color=cmap(14), label="TG")
        ax.plot(x_axis, TT, color=cmap(15), label="TT")

        ax.legend()

        plt.savefig(
            self.file_path.parent
            / Path(f"{self.file_path.name[:11]}_{motif_id}_dn_composition.png")
        )

        return

    def fit_probound(self):
        pass
