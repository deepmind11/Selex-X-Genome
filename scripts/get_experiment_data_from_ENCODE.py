import json
import random
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import requests


def save_experiment_data_as_json(search_result, base_dir):
    """
    Queries the ENCODE DB using the experiment accession and saves the json data at {target_dir}.

    Args:
        search_result (dict: @graph): An element of the @graph list returned by ENCODE search
        base_dir (Path):       The base directory where all the experiments will be saved.
    """

    expr_accession = search_result["accession"]

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
        [expr_data["possible_controls"][i]["accession"]]
        for i in range(len(expr_data["possible_controls"]))
    )

    with log_file.open(mode="r") as log:

        log.write(
            f"The query paramerters are tf={args.tf} organism={args.organism} expr_accession={expr_accession}{chr(10)}"
        )
        log.write(f"Date: {datetime.now():%c}{chr(10)}")
        log.write(f"List of Libraries: {*libraries,}{chr(10)}")
        log.write(f"Possible Controls: {*controls,}{chr(10)}")

    return

def download_and_process_fastq_files_for_experiment(experiment):
    


def save_control_data_for_experiment(experiment_json):
    """
    Takes as input an exepriment_data.json file. Check if the control experiment data exists, if not
    download it. And then proceed to download all the relevant files for that control experiment.

    Args:
        experiment_json (file(.json)): data/tf_organism/experiments/ENCSRXXXXXX/experiment_data.json
    """


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

        # Saving the experiment json.
        save_experiment_data_as_json(search_result, base_dir)

        # Starting a new process to download the files for that experiment.

        # !Downlaod the files for that experiment

        # !Downlaod the control files for that experiment
