#!/usr/bin/env python3

"""
Given a TF and organism. Downloads and transforms all the files from ENCODE.
The script is followed by the creation of the count tables.
"""


import subprocess
import sys
from pathlib import Path

sys.path.append("/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source")
from concurrent.futures import ThreadPoolExecutor

from base import DiskFile, ENCODE_Object, Experiment
from search import EncodeSearch

from source.diskfiles.fastq import PE_Fastq, SE_Fastq
from source.diskfiles.slurmjob import Slurmjob
from source.ENCODE.experiment import Control, TFChipSeq
from source.ENCODE.files import PE_File, SE_File
from source.ENCODE.library import Library

Human_MYC = EncodeSearch("CTCF", "Homo+sapiens", limit="all")
Human_MYC.search()

# print(Human_MYC.search_result)
Experiments = Human_MYC.get_experiments()


def process_controls(controls: list):
    """Downloads all the fastq files for the given set of control experiments. Transform the fastq files to fasta.
        Delete all fastq files.

    Args:
        controls (list): List of encode Control Objects.
    """

    for control in controls:
        try:
            control.fetchData()
            libraries = control.get_libraries()

            # Skip if control already exists
            # if Path("/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/Control/{control.accession}/").exists():
            #     continue

            for library in libraries:
                files = library.get_Files()

                for file in files:
                    # if fasta file exists then skip loop
                    if Path(
                        f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/Control/{control.accession}/{library.accession}/{file.accession}.fasta"
                    ).exists():
                        continue
                    fq = file.download(
                        Path(
                            f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/Control/{control.accession}/{library.accession}"
                        )
                    )
                    sub_fq = fq.subsample()
                    fq.delete()
                    sub_fq.slurm_transform()

        except Exception as e:
            print(f"Error: {control.accession}")
            print(e)


def process_exp(experiment: TFChipSeq):
    """For a given TFChipSeq object. It download all the important files along with the controls.
        It then process the fastq files. And calls a slurm job to transform it to fasta.
        Delete all fastq files.

    Args:
        experiment (TFChipSeq): A Chip seq ENCODE experiment.
    """
    try:
        experiment.fetchData()
        libraries = experiment.get_libraries()
        controls = experiment.get_controls()
        # Process the controls
        process_controls(controls)

        for library in libraries:
            files = library.get_Files()

            for file in files:
                # if fasta file exists then skip loop
                if Path(
                    f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/CTCF_Human/{experiment.accession}/{library.accession}/{file.accession}.fasta"
                ).exists():
                    continue

                fq = file.download(
                    Path(
                        f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/CTCF_Human/{experiment.accession}/{library.accession}"
                    )
                )
                sub_fq = fq.subsample()
                fq.delete()
                sub_fq.slurm_transform()

    except Exception as e:
        print(f"Error: {experiment.accession}")
        print(e)


# for experiment in Experiments:
#     process_exp(experiment)

if __name__ == "__main__":
    with ThreadPoolExecutor() as executor:
        executor.map(process_exp, Experiments)
