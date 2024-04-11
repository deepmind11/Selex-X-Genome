"""Gets all the fits from motif central"""

import json


def load_json(filename):
    with open(filename) as f:
        return json.load(f)


MOTIFCENTRAL = load_json(
    "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/motifcentral.json"
)
