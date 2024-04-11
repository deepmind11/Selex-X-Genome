#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path

sys.path.append("/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source")
import json
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from base import DiskFile, ENCODE_Object, Experiment
from search import EncodeSearch

from source.diskfiles.countTables import CountTable
from source.diskfiles.fastq import PE_Fastq, SE_Fastq
from source.diskfiles.slurmjob import Slurmjob
from source.ENCODE.experimentODE.experiment import Control, TFChipSeq
from source.ENCODE.files import PE_File, SE_File
from source.ENCODE.library import Library
from source.motifs.motif import Mononucleotide
from source.motifs.parse_motifcentral_json import MOTIFCENTRAL

max_motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[12])


Human_MYC = EncodeSearch("MYC", "Homo+sapiens", limit="all")
Human_MYC.search()

# print(Human_MYC.search_result)
Experiments = Human_MYC.get_experiments()


def get_fasta(control: Control):

    control_dr = Path(
        f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/Control/{control.accession}"
    )
    return list(control_dr.rglob("*.fasta"))


def create_table(experiment: TFChipSeq):

    try:
        experiment.fetchData()
        controls = experiment.get_controls()
        # For now let's select one control for now.
        control_fasta = get_fasta(controls[0])
        control_fasta = control_fasta[0]

        experiment_scores = {}

        libraries = experiment.get_libraries()

        for library in libraries:
            files = library.get_Files()

            for file in files:
                # if fasta file exists then skip loop
                r1_fasta = Path(
                    f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/MYC_Human/{experiment.accession}/{library.accession}/{file.accession}.fasta"
                )
                count_table_path = r1_fasta.parent / Path(f"{file.accession}.tsv")

                # Create the count table
                count_table = CountTable.create_from_fasta(
                    control_fasta, r1_fasta, count_table_path
                )
                count_table.processTable()
                count_table.zip()

                experiment_scores[(library.accession, file.accession)] = (
                    count_table.score(max_motif)
                )
                print(experiment.accession)
                print(experiment_scores[(library.accession, file.accession)])
        return experiment_scores

    except Exception as e:
        print(f"Error: {experiment.accession}")
        print(e)


if __name__ == "__main__":

    # Create a pool of processes
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        # Map the function to the list of items and collect the results
        results = pool.map(create_table, Experiments)

        # saving the results as json
        results_json = {}
        for i, experiment in enumerate(Experiments):
            results_json[experiment.accession] = results[i]

        json_file_path = Path(
            "/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/MYC_Human/scores.json"
        )
        json_file_path.touch(exist_ok=True)
        with open(json_file_path, "w") as file:
            json.dump(results_json, file, indent=4)
