"""
This program queries the ENCODE DB using its REST API. 
Saves all the data for an experiment in a json. data/tf_organism/ENCSRXXXXXX/experiment_data.json
Saves data for the control experiment as a json. data/control/ENCSRXXXXXX/exepriment_data.json
"""


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
