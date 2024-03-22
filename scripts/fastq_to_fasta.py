"""
Each TF+Organism pair has a bunch of experiments.
Every experiment has a bunch of libraries (SE or PE). These libraries contain fastq files. 
This script transfroms them into fasta taking SE/PE into account. Now, we are ready to build 
the count tables.
"""

from pathlib import Path
import json
import subprocess
from argparse import ArgumentParser

def slurm_job_to_transform_library(library, pe_flag=False):
    """
    Submits the {se/pe}_processing_template.sh scripts to the cluster.
    The above script takes the fastq files (full path) as input.
    
    It converts the fastq file to fasta. Finally, all intermediary files are deleted


    1. Create a ".job" file specifying the slurm configurations.
    2. Save the job file inside the library folder.
    3. Submit the job to the cluster.

    Args:
        library (Path): Path to library folder.
        pe_flag (bool): True if the library is paired-ended, False (default) otherwise.
    """
    
    lib_name = library.stem
    # Create a .job file
    library_transform_job_file = library / Path(
        f"{lib_name}_transform.job"
    )
    library_transform_job_file.touch()

    # Getting the fastq files
    if pe_flag:
        for file in library.iterdir():
            if file.suffix == ".fq" and '_1' in file.name:
                read1 = file
            elif file.suffix == ".fq" and '_2' in file.name:
                read2 = file
    else:
        for file in library.iterdir():
            if file.suffix == ".fq":
                read = file


    # Specify the configuration
    with library_transform_job_file.open(mode="w") as jf:
        jf.write("#!/bin/bash\n")
        jf.writelines(f"#SBATCH --job-name={lib_name}.job\n")
        jf.writelines(
            f"#SBATCH --output={str(library)}/{lib_name}_transform.out\n"
        )
        jf.writelines(
            f"#SBATCH --error={str(library)}/{lib_name}_transform.err\n"
        )
        jf.writelines("#SBATCH -c 1\n")
        jf.writelines("#SBATCH --mem-per-cpu=5G\n")
        jf.writelines("#SBATCH --account=hblab\n")
        jf.writelines("#SBATCH -t 2:30:00\n\n")
        if pe_flag:
            jf.writelines(
                f'{str(Path(__file__).parent/Path("pe_processing_template.sh"))} {str(read1)} {str(read2)}'
            )
        else:
            jf.writelines(
                f'{str(Path(__file__).parent/Path("se_processing_template.sh"))} {str(read)}'
            )
    # Submit Job to cluster
    subprocess.run(["sbatch", str(library_transform_job_file)])

    return


def fastq_to_fasta(experiment_json):
    """
    1. Iterates through each library in the experiment
    2. Converts the library into fasta.
    4. Looks for the control files for the experiment
    5. Checks if there is a fasta file
    6. If not transforms the control experiment to a fasta too. (Can be called recursively, as the control's controls is an empty list)

    Args:
        experiment_json (Path): Path to the experiment_data json file.
    """
    # Loading the json data in python dict.
    with open(experiment_json, "r") as f:
        expr_data = json.load(f)

    # Libraries of this experiment
    libraries = list(
        [
            expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(expr_data["replicates"]))
        ]
    )

    # Check if experiment is paired-end (set pe flag)
    pe_flag = (
        expr_data["files"][0]["run_type"]
        == "paired-ended"
    )
   
    # Transform the library
    for library in libraries:
        slurm_job_to_transform_library(experiment_json.parent / Path(library), pe_flag)


    # Controls for this experiment
    controls = list(
    [
        expr_data["possible_controls"][i]["accession"]
        for i in range(len(expr_data["possible_controls"]))
    ]
    )
    
    # Check if control fasta exists
    for control in controls:
        # Assume that fasta does not exist
        control_fasta_exist = False
        control_dr = Path(__file__).parent.parent / Path(
            f'data/Control/{control}')
        # Test that assumption
        for file in control_dr.glob("*"):
            if file.suffix == ".fasta":
                control_fasta_exist = True
                break
        # Check again if fasta file exists
        if not control_fasta_exist:
            fastq_to_fasta(control_dr /  Path("experiment_data.json"))
       

    return


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("tf", type=str)
    parser.add_argument("organism", type=str)

    args = parser.parse_args()

    experiment_DR = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}/experiments'
    )

    for experiment in experiment_DR.iterdir():
        if experiment.is_dir():
            experiment_json = experiment / Path("experiment_data.json")
            fastq_to_fasta(experiment_json)