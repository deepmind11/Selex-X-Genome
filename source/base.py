import requests


class ENCODE_Object:
    """Base class for different ENCODE objects (ENCFF, ENCSR, ENCLB, etc.)"""

    def __init__(self, accession):
        self.accession = accession

    def get_url(self):
        """Returns the url for the ENCODE object"""
        return f"https://www.encodeproject.org/{self.accession}"
    

    
