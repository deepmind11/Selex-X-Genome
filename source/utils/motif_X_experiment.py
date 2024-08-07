#!/usr/bin/env python3

import sqlite3
from argparse import ArgumentParser
from pathlib import Path

from diskfiles.countTables import CountTable
from ENCODE.experiment import TFChipSeq
from motifs.motif import Mononucleotide

if __name__ == "__main__":

    # Parsing the arguments
    parser = ArgumentParser()
    parser.add_argument("TF", type=str)
    parser.add_argument("Search_TF", type=str)
    parser.add_argument("organism", type=str)
    parser.add_argument("experiment", type=str)
    parser.add_argument("data_path", type=str)
    parser.add_argument("fit_id", type=str)
    parser.add_argument(
        "--psam",
        nargs="*",
        type=float,  # any type/callable can be used here
        default=[],
    )

    args = parser.parse_args()

    psam = args.psam
    data_path = args.data_path
    search_tf = args.Search_TF
    fit_id = args.fit_id
    motif = Mononucleotide(args.TF, args.organism, psam, fit_id)

    experiment = TFChipSeq(args.experiment)
    experiment.update_database()  # Update DB
    libraries = experiment.get_libraries()
    for library in libraries:
        library.update_database()  # Update DB
        file = library.get_Files()[0]  # By default returns a list of size 1
        file.update_database()  # Update DB
        # Check whether counttable exists
        if Path(
            f'{data_path}{search_tf}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}/{file.accession}.tsv.gz'
        ).exists():
            count_table = CountTable(
                Path(
                    f'{data_path}{search_tf}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}/{file.accession}.tsv.gz'
                )
            )
            count_table.update_database()
            count_table.score(motif, search_tf)
        else:
            fastq = file.download(
                Path(
                    f'{data_path}{search_tf}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}'
                )
            )
            subsampled_fq = fastq.subsample(size=1000000)
            fasta = subsampled_fq.transform_to_fasta()
            count_table = fasta.build_count_table(data_path)
            count_table.processTable()
            count_table.update_database()
            count_table.score(motif, search_tf)
