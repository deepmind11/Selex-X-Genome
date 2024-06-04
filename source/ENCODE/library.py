import sqlite3

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
        files: list[
            dict
        ],  # file dicts belonging to this library # ! All files must be fastq files associated with the library for
        # ! The methods to work!
        control: bool,
    ):
        super().__init__(accession)
        self.antibody = antibody
        self.technical_replicate_number = technical_replicate_number
        self.biological_replicate_number = biological_replicate_number
        self.experiment = experiment
        self.biosample = biosample
        self.run_type = run_type
        self.files = files
        self.control = control

    @classmethod
    def create_from_ENCODE_dict(
        cls, lib_dict: dict, run_type: RunType, lib_files: list[dict], control: bool
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
            control,
        )

    def get_Files(self, all: bool = False) -> list[SE_File] | list[PE_File]:
        """
        Returns a list of experiment associated fastq files as ENCODE File objects.
        If all is true returns all files, else returns the first one.
        """

        if all:
            count = len(self.files)
        else:
            count = 1

        if self.run_type == RunType.SE:
            return [
                SE_File.create_from_ENCODE_dict(
                    file,
                    RunType.SE,
                    self.experiment,
                    self.biosample,
                    self.technical_replicate_number,
                    self.biological_replicate_number,
                    self.control,
                    self.antibody,
                )
                for file in self.files
            ][:count]
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
                        SE_File.create_from_ENCODE_dict(
                            r1,
                            RunType.PE,
                            self.experiment,
                            self.biosample,
                            self.technical_replicate_number,
                            self.biological_replicate_number,
                            self.control,
                            self.antibody,
                        ),
                        SE_File.create_from_ENCODE_dict(
                            r2,
                            RunType.PE,
                            self.experiment,
                            self.biosample,
                            self.technical_replicate_number,
                            self.biological_replicate_number,
                            self.control,
                            self.antibody,
                        ),
                    )
                )

            return pe_files[:count]

    def update_database(self):
        "Adds information to the Selex_X_Genome.db database"
        # Connect to the database
        with sqlite3.connect(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            # Get the experiment_id from the experiments table
            try:
                cursor.execute(
                    """
                    SELECT id FROM experiments WHERE accession = ?
                    """,
                    (self.experiment,),
                )
                experiment_id = cursor.fetchone()[0]
            except:
                raise ValueError("Could not find experiment in Database")

            library = [
                (
                    self.accession,
                    self.antibody,
                    self.biosample,
                    self.technical_replicate_number,
                    self.biological_replicate_number,
                    experiment_id,
                )
            ]
            try:
                cursor.executemany(
                    """
                INSERT INTO "libraries" ("accession", "antibody", "biosample", "technical_rep_number", "biological_rep_number", "experiment_id") VALUES (?, ?, ?, ?, ?, ?)
                """,
                    library,
                )
                conn.commit()
            except sqlite3.IntegrityError:
                pass
            except Exception as e:
                print("Error:", e)
