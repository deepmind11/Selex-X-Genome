"""
This program queries the ENCODE DB using its REST API. 
Saves the search result as a json. data/encode_search/tf_organism.json
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


# ! An observation I made. So Bioreplicate * Technical Replicate -> Unique Library -> Ideally have one fastq files associated with it
# ! What does it mean if I find multiple such files? Ans) The reads come from multiple lanes that is why they are split up.
# ! 1) They could be paired ended (This info should be available in the exp metadata). 
# ! 2) The library is split into multiple fastq files.
# ! Maybe ask this on BioStars.


if __name__=="main":
    



def get_fastq_accession(search_object):
    """Returns the fastq accession for a given search_object (compenent of "@graph").
    Technically, a search object corresponds to an experiment accession.

    Args:
        search_object (dict: @graph): An element of the @graph list returned by ENCODE search
    Returns:
        fastq file accession for that experiment.
    """


test = ENCODE_search("CTCF", "Mus+musculus", "2")

print(json.dumps(test, indent=4))
