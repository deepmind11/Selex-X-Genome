from .base import ENCODE_Object
from .files import PE_File, SE_File


class Library(ENCODE_Object):
    """Class for ENCODE library object"""

    def __init__(
        self,
        accession,
        antibody,
        technical_replicate_number,
        biological_replicate_number,
        experiment,
        biosample,
        files: PE_File | SE_File,
    ):
        self.accession = accession
        self.antibody = antibody
        self.technical_replicate_number = technical_replicate_number
        self.biological_replicate_number = biological_replicate_number
        self.experiment = experiment
        self.biosample = biosample
        self.files = files

    # Constructor method from 'replicate' object and list of files for exp
    # Some attributes i need to set with the init method, but then other attributes I want to set by calling
    # For example self.files .

    def create_files(self) -> list[PE_File | SE_File]:
        return


# ('antibody', dict),
#  ('biological_replicate_number', int),
#  ('technical_replicate_number', int),
#  ('experiment', str),
#  ['library']['biosample'] # Info about the biosample.
#  list of files (files attribute.)
