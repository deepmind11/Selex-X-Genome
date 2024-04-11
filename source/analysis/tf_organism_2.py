#!/usr/bin/env python3

"""
Creates all the count tables, deletes the fasta files, scores the cnt table,
and returns enrichment of top bin. It saves a scores.json file at data/tf_organism. 
Which contains a sommary of scores for all the experiments.
"""

import json
import multiprocessing
from argparse import ArgumentParser
from pathlib import Path

from diskfiles.countTables import CountTable
from ENCODE.experiment import Control, TFChipSeq
from ENCODE.search import EncodeSearch
from motifs.motif import Mononucleotide
from motifs.parse_motifcentral_json import MOTIFCENTRAL


def get_fasta(control: Control) -> list[Path]:
    """Returns a list of paths to all the fasta files"""

    control_dr = Path(
        f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/Control/{control.accession}"
    )
    return list(control_dr.rglob("*.fasta"))


def create_table_and_score(experiment: TFChipSeq) -> dict:

    try:
        # Query Encode
        experiment.fetchData()

        # Get all the controls. This will be round 0 in the table.
        controls = experiment.get_controls()
        # ! For now selecting the first fasta file of the first control.
        control_fasta = get_fasta(controls[0])
        control_fasta = control_fasta[0]

        # Get all the libraries. These are the round 1s in the table.
        libraries = experiment.get_libraries()

        # Score for each library will be added to this dictionary
        experiment_scores = {}

        # Create a count table for each library and score it.
        for library in libraries:

            # Get all the files associated with the library.
            # ! Ideally it should just be one. I want to avoid repetition with split fastq files
            files = library.get_Files()

            for file in files:
                
                # Path to supposed fasta file
                r1_fasta = Path(
                    f"/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}/{file.accession}.fasta"
                )
                count_table_path = r1_fasta.parent / Path(f"{file.accession}.tsv")

                # Create the count table
                count_table = CountTable.create_from_fasta(
                    control_fasta, r1_fasta, count_table_path
                )
                count_table.processTable()
                count_table.zip()

                # ! I should make MOTIF a function parameter.
                experiment_scores[f"{library.accession}, {file.accession}"] = (
                    count_table.score(MOTIF)
                )
                print(experiment.accession)
                print(experiment_scores[f"{library.accession}, {file.accession}"])
        return experiment_scores

    except Exception as e:
        print(f"Error: {experiment.accession}")
        print(e)


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("TF", type=str)
    parser.add_argument("organism", type=str) 
    parser.add_argument("motifcentral_index", type=int)
    args = parser.parse_args()

    MOTIF = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[args.motifcentral_index])

    # Check that the organism has been specified in the right format.
    if "+" not in args.organism:
        raise ValueError("Organism must be in the format of Homo+sapiens")

    # Querying the encode DB
    query = EncodeSearch(args.TF, args.organism, limit="all")
    query.search()

    # All the experiments for the query
    Experiments = query.get_experiments()

    
    # Create a pool of processes
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        
        # Map the function to the list of items and collect the results
        results = pool.map(create_table_and_score, Experiments)

        # saving the results as json
        results_json = {}
        for i, experiment in enumerate(Experiments):
            results_json[experiment.accession] = results[i]

        json_file_path = Path(
            f'/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/{args.motifcentral_index}_scores.json'
        )
        json_file_path.touch(exist_ok=True)
        with open(json_file_path, "w") as file:
            json.dump(results_json, file, indent=4)
