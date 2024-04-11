import json


def load_json(filename):
    with open(filename) as f:
        return json.load(f)


MOTIFCENTRAL = load_json("../data/motifcentral.json")
