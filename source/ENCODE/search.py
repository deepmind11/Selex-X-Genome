import requests
from ENCODE.experiment import TFChipSeq


class EncodeSearch:
    """Class for searching for ENCODE TF ChIPseq experiments"""

    def __init__(self, tf, organism, limit="all", search_result: list[dict] = None):
        self.tf = tf
        self.organism = organism
        self.limit = limit
        self.search_result = search_result

    def search(self):
        # Force return from the server in JSON format
        headers = {"accept": "application/json"}

        # This searches the ENCODE database for the phrase "bone chip"
        url = (
            f"https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq"
            f"&target.label={self.tf}&replicates.library.biosample.donor.organism.scientific_name={self.organism}&status=released&files.run_type=single-ended&limit={self.limit}"
        )

        # GET the search result
        response = requests.get(url, headers=headers)

        # Extract the JSON response as a python dictionary
        self.search_result = response.json()
        self.search_result = self.search_result["@graph"]

        return self.search_result

    def get_experiments(self):
        """Returns List of TF ChipSeq Experiments"""
        if self.search_result is None:
            raise Exception("Fetch data first")
        else:
            experiments = list(
                [
                    self.search_result[i]["accession"]
                    for i in range(len(self.search_result))
                ]
            )
            return [TFChipSeq(experiment) for experiment in experiments]
