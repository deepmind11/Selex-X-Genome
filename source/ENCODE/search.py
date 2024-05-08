import requests
from ENCODE.experiment import TFChipSeq


class EncodeSearchError(Exception):
    pass


class EncodeSearch:
    """Class for searching for ENCODE TF ChIPseq experiments"""

    tax_id_dict = {
        7227: "Drosophila+melanogaster",
        9606: "Homo+sapiens",
        10090: "Mus+musculus",
        94885: "Pantherophis+guttatus",
    }

    def __init__(
        self,
        tf: str,
        organism: str,
        limit: str = "all",
    ):
        self.tf = tf
        self.organism = organism
        self.limit = limit
        try:
            self.search_result = self.search()
        except EncodeSearchError:
            self.search_result = None

    def search(self) -> dict:
        # Force return from the server in JSON format
        headers = {"accept": "application/json"}

        # This searches the ENCODE database for the phrase "bone chip"
        url = (
            f"https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq"
            f"&target.label={self.tf}&replicates.library.biosample.donor.organism.scientific_name={self.organism}&status=released&files.run_type=single-ended&limit={self.limit}"
        )

        # GET the search result
        response = requests.get(url, headers=headers)

        if response.status_code != 200:
            raise EncodeSearchError(f"Search error for {self.tf} and {self.organism}")

        # Extract the JSON response as a python dictionary
        self.search_result = response.json()
        self.search_result = self.search_result.get("@graph")

        if self.search_result is None:
            raise EncodeSearchError(f"No hits for {self.tf} and {self.organism}")

        return self.search_result

    def get_experiments(self):
        """Returns List of TF ChipSeq Experiments"""
        if self.search_result is None:
            raise EncodeSearchError
        else:
            experiments = list(
                [
                    self.search_result[i]["accession"]
                    for i in range(len(self.search_result))
                ]
            )
            return [TFChipSeq(experiment) for experiment in experiments]

    @classmethod
    def search_using_MOTIFCENTRAL_json(cls, motif: dict, limit: str = "all"):
        """Create ENCODE search object from an element of MOTIFCENTRAL dict"""
        # Get the organism
        tax_id = motif["metadata"]["factors"][0].get("tax_id")
        organism = EncodeSearch.tax_id_dict.get(tax_id)
        # Get the TF
        tf = motif["metadata"]["factors"][0].get("gene_symbol")
        tf = tf.upper()
        # Get the search object.
        return EncodeSearch(tf, organism, limit)
