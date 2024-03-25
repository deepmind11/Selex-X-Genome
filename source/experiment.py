from .base import ENCODE_Object
import requests


# Should I have am Experiment ABC and then have a class for each type of experiment, such as RNAseq, ChIPseq, etc

class Experiment(ENCODE_Object):
    """Class for ENCODE experiment object"""
    # How to do this?    
    self._search_result=None
    
    def __init__(self, accession):
        self.accession = accession


    def fetchData(self):
        """Gets the data from ENCODE"""
        # ! How to call super method to get URL?
        # Force return from the server in JSON format
        headers = {"accept": "application/json"}

        response = requests.get(url, headers=headers)

        # Extract the JSON response as a python dictionary
        self.search_results = response.json()
        
    def get_controls(self):
        pass

    def get_other_meta_data(self):
        pass

    def get_libraries(self):
        pass

    

        











class Control(Experiment):