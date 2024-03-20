#!/usr/bin/env python3

"""
Given a TF and organism. It creates count tables for all the libraries in each experiment.
The set of libraries is crossed with the set of controls, and each possible count table is stored 
in experiment/ENCSRXXXXXX/ENCLBXXXXXX/x<ENCSRXXXXXX (control accession)>.tsv 

Only works on cluster!
Requires Probound at ~/Probound/pb
"""

import json
import subprocess
from argparse import ArgumentParser
from pathlib import Path


def create_table(library, controls):
    """
    Creates all the count tables for the library with the given controls

    Args:
        library (Path): Path to the library fasta file
        controls (List[Paths]): List of Paths to control fasta files
    """
    for control in controls:
        count_table_path = library.parent / Path(f"x{control.parent.parent.name}.tsv")
        with open(count_table_path, "w") as output_file:
            subprocess.run(
                ["/burg/home/hg2604/ProBound/pb", "make-table", str(control), str(library)],
                stdout=output_file,
            )


def create_tables_for_exp(experiment):
    """
    Creates all the count tables for the experiment(ENCSRXXXXXX)

    Args:
        experiment (Path): Path to the experiment directory
    """

    experiment_json = experiment / Path("experiment_data.json")
    # Loading the json
    with experiment_json.open(mode="r") as f:
        expr_data = json.load(f)

    # Getting all the libraries for the experiment
    libraries = list(
        [
            expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(expr_data["replicates"]))
        ]
    )

    # Getting a list of library fasta files
    libraries_fasta = []
    for library in libraries:
        for file in (experiment / Path(library)).iterdir():
            if file.suffix == ".fa":
                libraries_fasta.append(file)

    # Getting all the controls for the experiment
    controls = list(
        [
            expr_data["possible_controls"][i]["accession"]
            for i in range(len(expr_data["possible_controls"]))
        ]
    )

    # Getting a list of control fasta files
    controls_fasta = []
    for control in controls:
        control_dir = Path(__file__).parent.parent / Path(f"data/Control/{control}")
        for path in sorted(control_dir.rglob("*")):
            if path.suffix == ".fa":
                controls_fasta.append(path)

    # Creating all the count tables
    for library in libraries_fasta:
        create_table(library, controls_fasta)


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("tf", type=str)
    parser.add_argument("organism", type=str)

    args = parser.parse_args()

    # Get the list of experiments.
    experiments = []
    for exp_folder in (
        Path(__file__).parent.parent
        / Path(f"data/{args.tf}_{args.organism.replace('+', '_')}")
        / Path("experiments")
    ).iterdir():
        experiments.append(exp_folder)

    # ! Ensure that this runs in parallel, won't subprocess achive that
    for experiment in experiments:
        create_tables_for_exp(experiment)
