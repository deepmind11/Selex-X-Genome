from ENCODE.base import Experiment
from ENCODE.library import Library


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

    @staticmethod
    def get_run_type(expr_data):
        for file in expr_data["files"]:
            if file["file_format"] == "fastq":
                return file["run_type"]

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
            run_type = TFChipSeq.get_run_type(self.expr_data)
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

        return {
            "life_stage_age": self.expr_data.get("life_stage_age"),
            "perturbed": self.expr_data.get("perturbed"),
            "lab": self.expr_data.get("lab").get("@id"),
            "biosample_class": self.expr_data.get("biosample_ontology").get(
                "classification"
            ),
            "developmental_slims": self.expr_data.get("biosample_ontology").get(
                "developmental_slims"
            ),
            "system_slims": self.expr_data.get("biosample_ontology").get(
                "central nervous system"
            ),
            "organ_slims": self.expr_data.get("biosample_ontology").get("organ_slims"),
            "cell_slims": self.expr_data.get("biosample_ontology").get("cell_slims"),
        }


class Control(Experiment):
    def __init__(self, accession, expr_data: dict = None):
        super().__init__(accession, expr_data)

    @staticmethod
    def get_run_type(expr_data):
        for file in expr_data["files"]:
            if file["file_format"] == "fastq":
                return file["run_type"]

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
            run_type = Control.get_run_type(self.expr_data)
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
