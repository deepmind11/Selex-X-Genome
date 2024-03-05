"""
This program queries the ENCODE DB using its REST API. 
Downloads the relevant fastq files for the experiments.
"""

import json

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

    return search_results["@graph"]  # This key has the search results.


test = ENCODE_search("CTCF", "Mus+musculus", "2")

print(json.dumps(test, indent=4))
