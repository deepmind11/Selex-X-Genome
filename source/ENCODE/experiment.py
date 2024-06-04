import sqlite3

from ENCODE.base import Experiment, FetchExperimentDataFailure, RunType
from ENCODE.library import Library


class RunTypeError(Exception):
    "Run-type of the experiment is not single-ended or paired-ended"
    pass


class MoreThanOneControl(Exception):
    "Experiment has more than one control"
    pass


class NoFastqFiles(Exception):
    "Experiment has no fastq files"
    pass


class NoControl(Exception):
    "Experiment has no controls"


class TFChipSeq(Experiment):
    """Class for TF ChipSeq experiment object"""

    def __init__(self, accession: str, expr_data: dict = None):
        super().__init__(accession, expr_data)
        self.fetchData()
        self.control = False
        if self.expr_data is None:
            raise FetchExperimentDataFailure

    def get_controls(self) -> "Control":  # Becasue control is defined later
        "Returns all the controls for the experiment, ideally one!"
        if len(self.expr_data["possible_controls"]) == 0:
            raise NoControl
        elif len(self.expr_data["possible_controls"]) > 1:
            raise MoreThanOneControl
        else:
            return Control(
                self.expr_data["possible_controls"][0]["accession"], self.accession
            )

    @staticmethod
    def get_run_type(expr_data: dict) -> RunType:
        "Returns the run-type of the experiment"
        for file in expr_data["files"]:
            if file["file_format"] == "fastq":
                if file["run_type"] == "single-ended":
                    return RunType.SE
                elif file["run_type"] == "paired-ended":
                    return RunType.PE
                else:
                    raise RunTypeError

        # If no fastq files found.
        raise NoFastqFiles

    def get_libraries(self) -> list[Library]:
        """Returns Libraries"""
        # Getting all the library accessions
        libraries = [
            self.expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(self.expr_data["replicates"]))
        ]

        # All the fastq files for each library
        library_fastq_files = dict(zip(libraries, [[] for i in range(len(libraries))]))

        for file in self.expr_data["files"]:
            if (
                "replicate" in file.keys() and file["file_format"] == "fastq"
            ):  # Should be fastq not fastq.gz
                library_fastq_files[file["replicate"]["library"][11:-1]].append(file)

        # Create libraries in a for loop
        run_type = TFChipSeq.get_run_type(self.expr_data)
        Libraries = []
        for i in range(len(libraries)):
            Libraries.append(
                Library.create_from_ENCODE_dict(
                    self.expr_data["replicates"][i],
                    run_type,
                    library_fastq_files[libraries[i]],
                    self.control,
                )
            )
        return Libraries

    def get_other_meta_data(self) -> dict:
        "Get important metadata as a subdict"
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
                "system_slims"
            ),
            "organ_slims": self.expr_data.get("biosample_ontology").get("organ_slims"),
            "cell_slims": self.expr_data.get("biosample_ontology").get("cell_slims"),
        }

    def update_database(self):
        "Adds information to the Selex_X_Genome.db database"
        # Connect to the database
        with sqlite3.connect(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            meta_data = self.get_other_meta_data()

            try:
                experiment = [
                    (
                        self.accession,
                        "n",
                        meta_data.get("life_stage_age"),
                        meta_data.get("perturbed"),
                        meta_data.get("lab"),
                        meta_data.get("biosample_class"),
                        ",".join(meta_data.get("developmental_slims", [])),
                        ",".join(meta_data.get("system_slims", [])),
                        ",".join(meta_data.get("organ_slims", [])),
                        ",".join(meta_data.get("cell_slims", [])),
                    )
                ]
                cursor.executemany(
                    """
                INSERT INTO "experiments" ("accession", "control", "life_stage_age", "perturbed", "lab", "biosample_class", "developmental_slims", "system_slims", "organ_slims", "cell_slims") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                    experiment,
                )
                conn.commit()
            except sqlite3.IntegrityError:
                pass
            except Exception as e:
                print("Error:", e)


class Control(Experiment):
    def __init__(self, accession: str, experiment: str, expr_data: dict = None):
        super().__init__(accession, expr_data)
        self.experiment = experiment
        self.fetchData()
        self.control = True
        if self.expr_data is None:
            raise FetchExperimentDataFailure

    @staticmethod
    def get_run_type(expr_data) -> RunType:
        "Returns the run-type of the experiment"
        for file in expr_data["files"]:
            if file["file_format"] == "fastq":
                if file["run_type"] == "single-ended":
                    return RunType.SE
                elif file["run_type"] == "paired-ended":
                    return RunType.PE
                else:
                    raise RunTypeError

        # If no fastq files found.
        raise NoFastqFiles

    def get_libraries(self) -> list[Library]:
        """Returns Libraries"""
        # Getting all the library accessions
        libraries = [
            self.expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(self.expr_data["replicates"]))
        ]

        # All the fastq files for each library
        library_fastq_files = dict(zip(libraries, [[] for i in range(len(libraries))]))

        for file in self.expr_data["files"]:
            if (
                "replicate" in file.keys() and file["file_format"] == "fastq"
            ):  # Should be fastq not fastq.gz
                library_fastq_files[file["replicate"]["library"][11:-1]].append(file)

        # Create libraries in a for loop
        run_type = Control.get_run_type(self.expr_data)
        Libraries = []
        for i in range(len(libraries)):
            Libraries.append(
                Library.create_from_ENCODE_dict(
                    self.expr_data["replicates"][i],
                    run_type,
                    library_fastq_files[libraries[i]],
                    self.control,
                )
            )
        return Libraries

    def update_database(self):
        "Adds information to the Selex_X_Genome.db database"
        # Connect to the database
        with sqlite3.connect(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            try:
                experiment = [(self.accession, "y")]
                cursor.executemany(
                    """
                INSERT INTO "experiments" ("accession", "control") VALUES (?, ?)
                """,
                    experiment,
                )
                conn.commit()
            except:
                pass
