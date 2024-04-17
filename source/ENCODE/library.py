from ENCODE.base import ENCODE_Object, RunType
from ENCODE.files import PE_File, SE_File


class Library(ENCODE_Object):
    """Class for ENCODE library object"""

    def __init__(
        self,
        accession: str,
        antibody: str,
        technical_replicate_number: int,
        biological_replicate_number: int,
        experiment: str,
        biosample: str,
        run_type: RunType,
        files: list[dict],  # file dicts belonging to this library
    ):
        super().__init__(accession)
        self.antibody = antibody
        self.technical_replicate_number = technical_replicate_number
        self.biological_replicate_number = biological_replicate_number
        self.experiment = experiment
        self.biosample = biosample
        self.run_type = run_type
        self.files = files

    # Constructor method from 'replicate' object and list of files for exp
    # Some attributes i need to set with the init method, but then other attributes I want to set by calling
    # For example self.files .

    @classmethod
    def create_from_ENCODE_dict(
        cls, lib_dict: dict, run_type: RunType, lib_files: list[dict]
    ):
        """Create a Library object from a dictionary"""
        return Library(
            lib_dict.get("library", {}).get("accession"),
            lib_dict.get("antibody", {}).get("accession"),
            lib_dict.get("technical_replicate_number"),
            lib_dict.get("biological_replicate_number"),
            lib_dict.get("experiment")[13:-1],
            lib_dict.get("library", {}).get("biosample", {}).get("accession"),
            run_type,
            lib_files,
        )

    def get_Files(self) -> list[SE_File] | list[PE_File]:

        if self.run_type == RunType.SE:
            return [SE_File.create_from_ENCODE_dict(file) for file in self.files]
        elif self.run_type == RunType.PE:
            # Get a list of tuples, each containing a pair
            # ! Can this be done more elegantly
            files = self.files.copy()
            pe_files = []
            while len(files) != 0:
                r1 = files[0]
                files.remove(r1)
                for file in files:
                    if file["accession"] == r1["paired_with"]:
                        r2 = file
                        files.remove(file)
                        break
                pe_files.append(
                    PE_File(
                        SE_File.create_from_ENCODE_dict(r1),
                        SE_File.create_from_ENCODE_dict(r2),
                    )
                )

            return pe_files
