#!/usr/bin/env python3

import os
from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
from diskfiles.countTables import CountTable
from motifs.motif import Mononucleotide
from motifs.parse_motifcentral_json import MOTIFCENTRAL


def find(name, path="/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data"):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


if __name__ == "__main__":

    # Parsing the arguments
    parser = ArgumentParser()
    parser.add_argument("motif_number", type=str)
    parser.add_argument("file_acc", type=str)

    args = parser.parse_args()

    motif_number = int(args.motif_number)

    motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[motif_number])

    # searching for the file
    count_table_path = find(args.file_acc + ".tsv.gz")

    print(count_table_path)

    table = CountTable(Path(count_table_path))

    fig, ax = table.plot_enrichment_vs_bin(motif)

    plt.show()

    # To save the plot to a file
    fig.savefig(
        f"/burg/home/hg2604/hblab/Projects/Selex-X-Genome/enrichment_vs_bin/{args.file_acc}.png"
    )
