"""Classes for different types of motifs. Mono, Di, Coop, etc."""

import numpy as np
from motifs.base import Motif


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
