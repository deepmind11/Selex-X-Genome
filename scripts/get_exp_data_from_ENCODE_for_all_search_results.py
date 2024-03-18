"""
Iterates through all the search results for a given (TF + organism) from the 
get_search_results_from_ENCODE.py script.

For each hit, it queries the ENCODE DB and download the experiment json. /data/Tf+organism/experiment/ENCSRXXXXXX/experiment.json
It then downloads fastq files for all libraries in the experiment. Processes these fastq files
Then it download all the controls (if necessary) for the experiment. /data/Control/ENCSRXXXXXX

"""

import json
import subprocess
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import requests


def save_experiment_data_as_json(expr_accession, base_dir):
    """
    Queries the ENCODE DB using the experiment accession and saves the json data at {target_dir}.

    Args:
        expr_accession (String): The experiment accession ENCSRXXXXXX
        base_dir (Path):       The base directory where all the experiments will be saved.
    """

    # Force return from the server in JSON format
    headers = {"accept": "application/json"}

    # This searches the ENCODE database for the phrase "bone chip"
    url = f"https://www.encodeproject.org/experiments/{expr_accession}"

    # GET the search result
    response = requests.get(url, headers=headers)

    # Extract the JSON response as a python dictionary
    expr_data = response.json()

    # Constructing the path to save the json
    target_dir = base_dir / Path(expr_accession)
    target_dir.mkdir(parents=True, exist_ok=True)

    # Saving search results as json
    json_file = target_dir / Path("experiment_data.json")
    json_file.touch()

    with json_file.open(mode="w") as js:
        json.dump(expr_data, js, indent=4)

    # Logging Info. Storing all the libraries and the control belonging to an experiment
    log_file = target_dir / Path("experiment.log")
    log_file.touch()

    libraries = list(
        [
            expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(expr_data["replicates"]))
        ]
    )
    controls = list(
        [
            expr_data["possible_controls"][i]["accession"]
            for i in range(len(expr_data["possible_controls"]))
        ]
    )

    with log_file.open(mode="w") as log:

        log.write(
            f"The query paramerters are tf={args.tf} organism={args.organism} expr_accession={expr_accession}{chr(10)}"
        )
        log.write(f"Date: {datetime.now():%c}{chr(10)}")
        log.write(f"List of Libraries: {*libraries,}{chr(10)}")
        log.write(f"Possible Controls: {*controls,}{chr(10)}")

    return


def slurm_job_to_download_and_process_files(experiment_json):
    """
    Submits the download_and_subsample_fastq_files.py script to the cluster.
    The above script takes the experiment_json(Path) as an argument.
    It downloads the fastq_files for all the libraries, subsamples (def = 2 million),
    and converts the fastq file to fasta. Finally, all fastq file are deleted/


    1. Create a ".job" file specifying the slurm configurations.
    2. Save the job file to experiment_json.parent
    3. Submit the job to the cluster.

    Args:
        experiment_json (Path): Path to experiment json file.
    """
    # Getting the experiment accession
    with experiment_json.open(mode="r") as f:
        expr_data = json.load(f)
    expr_accession = expr_data["accession"]

    # Create a .job file
    download_and_process_job_file = experiment_json.parent / Path(
        f"{expr_accession}.job"
    )
    download_and_process_job_file.touch()

    # Specify the configuration
    with download_and_process_job_file.open(mode="w") as jf:
        jf.write("#!/bin/bash\n")
        jf.writelines(f"#SBATCH --job-name={expr_accession}.job\n")
        jf.writelines(
            f"#SBATCH --output={str(experiment_json.parent)}/{expr_accession}.out\n"
        )
        jf.writelines(
            f"#SBATCH --error={str(experiment_json.parent)}/{expr_accession}.err\n"
        )
        jf.writelines("#SBATCH -c 1\n")
        jf.writelines("#SBATCH --mem-per-cpu=5G\n")
        jf.writelines("#SBATCH --account=hblab\n")
        jf.writelines("#SBATCH -t 2:30:00\n\n")
        jf.writelines(
            f'{str(Path(__file__).parent/Path("download_and_subsample_fastq_files.py"))} {str(experiment_json)}'
        )
    # Submit Job to cluster
    subprocess.run(["sbatch", str(download_and_process_job_file)])

    return


def save_control_data_for_experiment(experiment_json):
    """
    Takes as input an exepriment_data.json file. Check if the control experiment data exists, if not
    download it. And then proceed to download all the relevant files for that control experiment.

    Args:
        experiment_json (file(.json)): data/tf_organism/experiments/ENCSRXXXXXX/experiment_data.json
    """
    # Loading the json
    with experiment_json.open(mode="r") as f:
        expr_data = json.load(f)

    # Getting the controls
    controls = list(
        [
            expr_data["possible_controls"][i]["accession"]
            for i in range(len(expr_data["possible_controls"]))
        ]
    )

    for control in controls:
        # Get the DR for the control experiment
        control_DR = Path(__file__).parent.parent / Path(f"data/Control/{control}")

        # Check if the control experiment exists; if it does then go to next control
        if control_DR.exists():
            continue
        else:
            control_DR.mkdir(parents=True)

            # Saving the JSON for the control experiment
            save_experiment_data_as_json(control, control_DR.parent)

            # Path to the control experiment json
            control_experiment_json = control_DR / Path("experiment_data.json")
            # Submitting job to cluster for downloading all the files
            slurm_job_to_download_and_process_files(control_experiment_json)

    return


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("tf", type=str)
    parser.add_argument("organism", type=str)

    args = parser.parse_args()

    # loading the search results
    search_result_json = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}/search_results.json'
    )

    with search_result_json.open(mode="r") as srj:
        search_results = json.load(srj)

    search_results = search_results["@graph"]

    # Querying the ENCODE DB for all the hits(expriments) and saving as json.
    base_dir = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}/experiments'
    )

    for search_result in search_results:

        # Getting the experiment accession
        expr_accession = search_result["accession"]

        # Saving the experiment json.
        save_experiment_data_as_json(expr_accession, base_dir)

        # Starting a new process to download the files for current experiment.
        # Getting the path to experiment json
        target_dir = base_dir / Path(expr_accession)
        experiment_json = target_dir / Path("experiment_data.json")
        # Submitting job to cluster
        slurm_job_to_download_and_process_files(experiment_json)

        # Downlaod the control files for that experiment
        save_control_data_for_experiment(experiment_json)
