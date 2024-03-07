"""
This program queries the ENCODE DB using its REST API. Takes as input the search results.

1. Saves all the data for an experiment in a json. data/tf_organism/experiments/ENCSRXXXXXX/experiment_data.json
2. Saves data for the control experiment as a json. data/control/ENCSRXXXXXX/exepriment_data.json

3. Saves only one (picks one randomly, if multiple available. Downloads 1 pair for PE) fastq files belong to a particular 
   library (defined by tech replicate # & bio replicate #) at
   data/tf_organism/experiments/ENCSRXXXXXX/ENCLBXXXXXX/file.fastq

   3.1) Subsamples the fastq file 2,000,000. gzips it
   3.2) deletes the downloaded file.
"""

import json
import logging
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import requests


def get_experiment_data_from_ENCODE(search_result):
    """
    Queries the ENCODE DB using the experiment accession and stores the json data in a python dict.

    Args:
        search_result (dict: @graph): An element of the @graph list returned by ENCODE search
    Returns:
        expr_data (dict): Data associated with a given experiment.
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

    return expr_data


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

    # Querying the ENCODE DB for all the hits(expriments)
    base_dir = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}/experiments'
    )

    for search_result in search_results:

        expr_data = get_experiment_data_from_ENCODE(search_result)
        expr_accession = search_result["accession"]

        target_dir = base_dir / Path(expr_accession)
        target_dir.mkdir(parents=True, exist_ok=True)

        # Saving search results as json
        json_file = target_dir / Path("experiment_data.json")
        json_file.touch()

        with json_file.open(mode="w") as js:
            json.dump(expr_data, js, indent=4)

        # Logging Info. Storing all the libraries and the control belonging to an experiment
        # Logging Info
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

        logging.basicConfig(
            filename=log_file,
            filemode="w",
            level=logging.DEBUG,
            format="%(name)s - %(levelname)s - %(message)s",
        )
        logging.info(
            f"The query paramerters are tf={args.tf} organism={args.organism} expr_accession={expr_accession}"
        )
        logging.info(f"Date: {datetime.now():%c}")
        logging.info(f"List of Libraries: {*libraries,}")
        logging.info(f"Possible Controls: {*controls,}")
