#!/usr/bin/env python3

import sqlite3
from argparse import ArgumentParser
from pathlib import Path

from ENCODE.experiment import TFChipSeq
from motifs.motif import Mononucleotide

if __name__ == "__main__":

    # Parsing the arguments
    parser = ArgumentParser()
    parser.add_argument("TF", type=str)
    parser.add_argument("organism", type=str)
    parser.add_argument("experiment", type=str)
    parser.add_argument("psam", type=str)  # ! Mononucleotide PSAM
    args = parser.parse_args()

    psam = args.psam.split(",")
    psam = list([int(x) for x in psam])

    motif = Mononucleotide(args.TF, args.organism, psam)

    experiment = TFChipSeq(args.experiment)
    libraries = experiment.get_libraries()
    for library in libraries:
        file = library.get_Files()[0]

        # Check whether counttable exists
        if Path(
            f'/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}/{file.accession}.tsv.gz'
        ).exists():
            count_table.score(motif)
        else:
            fastq = file.download(
                Path(
                    f'/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}'
                )
            )
            subsampled_fq = fastq.subsample(size=1000000)
            fasta = subsampled_fq.transform_to_fasta()
            count_table = fasta.build_count_table()
            count_table.processTable()
            scores = count_table.score(motif)

            # Add Score to DB
            with sqlite3.connect(
                "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
            ) as conn:
                cursor = conn.cursor()

                # Get the r1_file_id from the files table
                try:
                    cursor.execute(
                        """
                        SELECT "id" FROM "count_tables" WHERE "r1_file" =?
                        """,
                        (count_table.file_path.name[:11]),
                    )
                    count_table_id = cursor.fetchone()[0]
                except:
                    raise ValueError("Could not find count_table in Database")

                # Value to be added to DB
                motif_data = [
                    (
                        "Mono",
                        motif.tf,
                        motif.organism,
                        count_table_id,
                        scores.get("avg_score"),
                        scores.get("r0_count"),
                        scores.get("r1_count"),
                        scores.get("enrichment"),
                    )
                ]
                # Update DB
                try:
                    cursor.executemany(
                        """
                    INSERT INTO "count_tables" ("type", "tf", "organism", "count_table_id", "score", "r0_count", "r1_count", "enrichment") VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                        motif_data,
                    )
                    conn.commit()
                except:
                    pass
