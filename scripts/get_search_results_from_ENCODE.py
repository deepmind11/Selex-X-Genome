"""
This program queries the ENCODE DB using its REST API. Query = <GeneName>  + <Organism:Homo+Sapiens>
Saves the search result as a json. data/tf_organism/search_result.json
"""

import json
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

import requests


def ENCODE_search(tf, organism, limit="all"):
    """Searches the encode database for ChIPseq experiments for given tf and organism.

    Args:
        tf (string): The Gene Symbol for the TF. Ex: CTCF, EBF1, etc.
        organism (string): The scientific name of the organism. Spaces are +. Ex: Homo+Sapiens, Mus+muculus,etc
        limit (string): The number of results to return. Default is 'all'
    """

    # Force return from the server in JSON format
    headers = {"accept": "application/json"}

    # This searches the ENCODE database for the phrase "bone chip"
    url = (
        f"https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq"
        f"&target.label={tf}&replicates.library.biosample.donor.organism.scientific_name={organism}&status=released&limit={limit}"
    )

    # GET the search result
    response = requests.get(url, headers=headers)

    # Extract the JSON response as a python dictionary
    search_results = response.json()

    return search_results  # This key has the search results.


# ! An observation I made. So Bioreplicate * Technical Replicate -> Unique Library -> Ideally have one fastq files associated with it
# ! What does it mean if I find multiple such files? Ans) The reads come from multiple lanes that is why they are split up.
# ! 1) They could be paired ended (This info should be available in the exp metadata).
# ! 2) The library is split into multiple fastq files.
# ! Ask this on BioStars.


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("tf", type=str)
    parser.add_argument("organism", type=str)
    parser.add_argument("--limit", type=str, default="all")

    args = parser.parse_args()

    target_dir = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}'
    )
    target_dir.mkdir(parents=True, exist_ok=True)

    # Saving search results as json
    json_file = target_dir / Path("search_results.json")
    json_file.touch()
    search_results = ENCODE_search(args.tf, args.organism, args.limit)

    with json_file.open(mode="w") as js:
        json.dump(search_results, js, indent=4)

    # Saving meta data about search results into a log file
    log_file = target_dir / Path("search.log")
    log_file.touch()
    total_hits = len(search_results["@graph"])

    with log_file.open(mode="w") as f:
        f.write(
            f"The query paramerters are tf={args.tf} organism={args.organism} limit={args.limit}{chr(10)}"  # chr(10) => new line
        )
        f.write(f"Date: {datetime.now():%c}{chr(10)}")
        f.write(f"Total hits: {total_hits}{chr(10)}")
