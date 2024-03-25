import sys

import requests

sys.path.append("/Users/hgz/Documents/researchHBLab/Projects/Selex-X-Genome/source")
from base import ENCODE_Object, Experiment
from library import Library


class TFChipSeq(Experiment):
    """Class for TF ChipSeq experiment object"""

    def __init__(self, accession, expr_data: dict = None):
        super().__init__(accession, expr_data)

    def get_controls(self) -> "list[Control]":  # Becasue control is defined later
        controls = list(
            [
                self.expr_data["possible_controls"][i]["accession"]
                for i in range(len(self.expr_data["possible_controls"]))
            ]
        )
        return [Control(control) for control in controls]

    def get_libraries(self) -> list[Library]:
        """Returns Libraries"""
        if self.expr_data is None:
            raise Exception("Fetch data first")
        else:
            libraries = list(
                [
                    self.expr_data["replicates"][i]["library"]["accession"]
                    for i in range(len(self.expr_data["replicates"]))
                ]
            )

            library_fastq_files = dict(
                zip(libraries, [[] for i in range(len(libraries))])
            )

            for file in self.expr_data["files"]:
                if (
                    "replicate" in file.keys() and file["file_format"] == "fastq"
                ):  # Should be fastq not fastq.gz
                    library_fastq_files[file["replicate"]["library"][11:-1]].append(
                        file
                    )

            self.expr_data["replicates"]

            # Create libraries in a for loop
            run_type = self.expr_data["files"][0]["run_type"]
            Libraries = []
            for i in range(len(libraries)):
                Libraries.append(
                    Library.create_from_ENCODE_dict(
                        self.expr_data["replicates"][i],
                        run_type,
                        library_fastq_files[libraries[i]],
                    )
                )
            return Libraries

    def get_other_meta_data(self):
        pass


class Control(Experiment):
    def __init__(self, accession, experiment_json: dict = None):
        super().__init__(accession, experiment_json)

    def get_libraries(self) -> list[Library]:
        """Returns Libraries"""
        if self.expr_data is None:
            raise Exception("Fetch data first")
        else:
            libraries = list(
                [
                    self.expr_data["replicates"][i]["library"]["accession"]
                    for i in range(len(self.expr_data["replicates"]))
                ]
            )

            library_fastq_files = dict(
                zip(libraries, [[] for i in range(len(libraries))])
            )

            for file in self.expr_data["files"]:
                if (
                    "replicate" in file.keys() and file["file_format"] == "fastq"
                ):  # Should be fastq not fastq.gz
                    library_fastq_files[file["replicate"]["library"][11:-1]].append(
                        file
                    )

            self.expr_data["replicates"]

            # Create libraries in a for loop
            run_type = self.expr_data["files"][0]["run_type"]
            Libraries = []
            for i in range(len(libraries)):
                Libraries.append(
                    Library.create_from_ENCODE_dict(
                        self.expr_data["replicates"][i],
                        run_type,
                        library_fastq_files[libraries[i]],
                    )
                )
            return Libraries
