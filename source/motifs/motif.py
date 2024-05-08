"""Classes for different types of motifs. Mono, Di, Coop, etc."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ENCODE.search import EncodeSearch

import numpy as np
import pandas as pd
from motifs.base import Motif
from motifs.parse_motifcentral_json import MOTIFCENTRAL
from utils.slurmjob import Slurmjob

# Add methods to get all the ENCODE files for a particular TF (not necessarily same) (Use the ENCODESearch Object)
# Build count-tables (parametrised by the (-100,200) cuts) (Don't build if already exist)
# Score the count tables
# Add data to SQL DB


class Mononucleotide(Motif):
    def __init__(self, tf, organism, psam: list):
        super().__init__(tf, organism)
        self.psam = psam

    @classmethod
    def create_from_motif_central(cls, motifcentral: dict):  # ! Make this more robust
        """Instantiate mononucleotide motif from motifcentral model json"""
        tf = motifcentral["metadata"]["factors"][0]["gene_symbol"]
        organism = motifcentral["metadata"]["factors"][0]["tax_id"]
        psam = motifcentral["coefficients"]["bindingModes"][0]["mononucleotide"]
        return cls(tf, organism, psam)

    @staticmethod
    def score_window(numpy_motif, seq):
        """Scores a window the size of motif"""

        score = 0
        for i, base in enumerate(seq):
            match base:
                case "A":
                    score += numpy_motif[0, i]
                case "C":
                    score += numpy_motif[1, i]
                case "G":
                    score += numpy_motif[2, i]
                case "T":
                    score += numpy_motif[3, i]

        return score

    def score_seq(self, seq: str):
        # Scores a sequence. Assert sequence lenght > motif lenght
        motif_size = int(len(self.psam) / 4)
        numpy_motif = np.array(self.psam).reshape(4, motif_size, order="F")

        score = 0
        for i in range(0, len(seq) - 2):
            score += Mononucleotide.score_window(numpy_motif, seq[i : i + motif_size])

        return score

    def score_random_kmers(
        self,
        data: Path = Path(
            "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/random_sequences.csv"
        ),
    ):
        """Get the average score and score distribution of the motif scored against data/random_sequences.csv"""
        df = pd.read_csv(data, header=None, index_col=False)
        df.columns = ["seq"]
        df["score"] = df["seq"].apply(lambda x: self.score_seq(x))
        # print(df["score"].mean())
        print(df)
        ax = df["score"].hist()
        # s is an instance of Series
        fig = ax.get_figure()
        fig.savefig("/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/figure.png")

        return

    def motif_X_TF(self, searchResult: EncodeSearch):

        experiments = searchResult.get_experiments()

        PSAM = list([str(i) for i in self.psam])

        for experiment in experiments:

            # Create Slurm Job
            job_path = Path(
                f'/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/{searchResult.tf}_{searchResult.organism.replace("+","_")}/{experiment.accession}/{searchResult.tf}_Score.job'
            )
            job_path
            slurmjob = Slurmjob(
                file_path=job_path,
                job_name=f'Score_{searchResult.tf}_{searchResult.organism.replace("+","_")}',
                job_script=Path(
                    "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/source/utils/motif_X_experiment.py"
                ),
                job_params=(
                    searchResult.tf,
                    searchResult.organism,
                    experiment.accession,
                    ",".join(PSAM),
                ),
                output=f'Score_{searchResult.tf}_{searchResult.organism.replace("+","_")}',
                cores=4,
                time=8,
            )
            slurmjob.submitJob()
