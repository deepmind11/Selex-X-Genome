#!/usr/bin/env python3

"""This has already been specifically implemented for myc and ctcf in "human_ctcf_cnt_tables.py" """

import subprocess
import sys
from argparse import ArgumentParser
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

max_motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[222])


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
                    f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/CTCF_Human/{experiment.accession}/{library.accession}/{file.accession}.fasta"
                )
                count_table_path = r1_fasta.parent / Path(f"{file.accession}.tsv")

                # Create the count table
                count_table = CountTable.create_from_fasta(
                    control_fasta, r1_fasta, count_table_path
                )
                count_table.processTable()
                count_table.zip()

                experiment_scores[f"{library.accession}, {file.accession}"] = (
                    count_table.score(max_motif)
                )
                print(experiment.accession)
                print(experiment_scores[f"{library.accession}, {file.accession}"])
        return experiment_scores

    except Exception as e:
        print(f"Error: {experiment.accession}")
        print(e)


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("experiment", type=str)

    args = parser.parse_args()

    expr_acc = args.experiment

    expr_obj = TFChipSeq(expr_acc)
    expr_obj.fetchData()

    result = create_table(expr_obj)

    result_path = Path(
        f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/CTCF_Human/{expr_acc}/experiment_scores.json"
    )
    result_path.touch(exist_ok=True)
    with open(result_path, "w") as file:
        json.dump(result, file, indent=4)
