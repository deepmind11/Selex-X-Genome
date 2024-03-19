#!/usr/bin/env python3

"""
Given a TF and organism. It creates count tables for all the libraries in each experiment.
The set of libraries is crossed with the set of controls, and each possible count table is stored 
in experiment/ENCSRXXXXXX/ENCLBXXXXXX/x<ENCSRXXXXXX (control accession)>.tsv 

Only works on cluster!
Requires Probound at ~/Probound/pb
"""

import subprocess
from pathlib import Path


def create_table(library, controls):
    """
    Creates all the count tables for the library with the given controls

    Args:
        library (Path): Path to the library fasta file
        controls (List[Paths]): List of Paths to control fasta files
    """
    for control in controls:
        count_table_path = library.parent / Path(f"x{control.stem}.tsv")
        subprocess.run(
            ["~/Probound/pb", "make-table", str(control), str(library)],
            ">",
            str(count_table_path),
        )


def create_tables_for_experiment(experiment):
    """
    Creates all the count tables for the experiment(ENCSRXXXXXX)

    Args:
        experiment (_type_): _description_
    """
