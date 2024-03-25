import requests


class EncodeSearch:
    """Class for searching for ENCODE experiments"""

    def __init__(self, tf, organism, limit="all"):
        self.tf = tf
        self.organism = organism
        self.limit = limit

    def search(self):
        # Force return from the server in JSON format
        headers = {"accept": "application/json"}

        # This searches the ENCODE database for the phrase "bone chip"
        url = (
            f"https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq"
            f"&target.label={self.tf}&replicates.library.biosample.donor.organism.scientific_name={self.organism}&status=released&limit={self.limit}"
        )

        # GET the search result
        response = requests.get(url, headers=headers)

        # Extract the JSON response as a python dictionary
        search_results = response.json()

        return search_results["@graph"]
