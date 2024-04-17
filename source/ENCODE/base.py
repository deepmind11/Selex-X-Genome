from enum import Enum

import requests


class FetchExperimentDataFailure(Exception):
    "Error to be raised when Fetch Data gets a 400 response code"
    pass


class RunType(Enum):
    """Enum for the different run types"""

    SE = "single-ended"
    PE = "paired-ended"


class ENCODE_Object:
    """Base class for different ENCODE objects (ENCFF, ENCSR, ENCLB, etc.)"""

    def __init__(self, accession: str):
        self.accession = accession

    def get_url(self) -> str:
        """Returns the url for the ENCODE object"""
        return f"https://www.encodeproject.org/{self.accession}"


class Experiment(ENCODE_Object):
    """Class for ENCODE experiment object"""

    def __init__(self, accession: str, expr_data: dict = None):
        super().__init__(accession)
        self.expr_data = expr_data

    def fetchData(self) -> dict:
        """Gets the data from ENCODE"""
        if self.expr_data is not None:
            return self.expr_data
        else:
            # Get the URL
            url = self.get_url()
            # Force return from the server in JSON format
            headers = {"accept": "application/json"}
            # Query the DB
            response = requests.get(url, headers=headers)
            if response.status_code != 200:
                raise FetchExperimentDataFailure
            else:
                # Extract the JSON response as a python dictionary
                self.expr_data = response.json()
                return self.expr_data
