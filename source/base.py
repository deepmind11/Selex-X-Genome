from pathlib import Path

import requests


class ENCODE_Object:
    """Base class for different ENCODE objects (ENCFF, ENCSR, ENCLB, etc.)"""

    def __init__(self, accession):
        self.accession = accession

    def get_url(self):
        """Returns the url for the ENCODE object"""
        return f"https://www.encodeproject.org/{self.accession}"


class DiskFile:
    """Base class for files on disk"""

    def __init__(
        self, file_path: Path = None
    ):  # Polymorphism in the tpyes supported by file_path. For ex: Path | list(Path)
        self.file_path = file_path

    def delete(self):
        if self.file_path is not None:
            self.file_path.unlink()
            self.file_path = None


class Experiment(ENCODE_Object):
    """Class for ENCODE experiment object"""

    def __init__(self, accession, expr_data: dict = None):
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
            # Extract the JSON response as a python dictionary
            self.expr_data = response.json()

            return self.expr_data
