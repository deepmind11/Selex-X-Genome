from __future__ import annotations

import sqlite3
from pathlib import Path

import requests
from diskfiles.fastq import PE_Fastq, SE_Fastq
from ENCODE.base import ENCODE_Object, RunType


class NoFileAvailable(Exception):
    """Exception for when no file is available"""

    pass


class FileDownloadError(Exception):
    pass


class SE_File(ENCODE_Object):
    """Class for a SE ENCODE file (ENCFFXXXXXX) object"""

    def __init__(
        self,
        accession: str,
        read_count: int,
        file_format: str,
        no_file_available: bool,
        platform: str,
        read_length: int,
        experiment: str,
        control: bool,
        library: str,  # Belongs to
        biosample: str,  # Sample from which reads were derived
        technical_replicate_number: int,
        biological_replicate_number: int,
        antibody: str,
        href: str,
        run_type: RunType,  # This has two optional parameters pe or se
        paired_end=None,  # Will be none for SE (Assert this)
        paired_with: str = None,  # WILL be none for SE (Assert this)
    ):
        super().__init__(accession)
        self.read_count = read_count
        self.file_format = file_format
        self.no_file_available = no_file_available
        self.platform = platform
        self.read_length = read_length
        self.experiment = experiment
        self.control = control
        self.library = library
        self.biosample = biosample
        self.technical_replicate_number = technical_replicate_number
        self.biological_replicate_number = biological_replicate_number
        self.antibody = antibody
        self.href = href
        self.run_type = run_type
        self.paired_end = paired_end
        self.paired_with = paired_with

        # Assert
        if self.run_type == RunType.SE:
            assert self.paired_end is None and self.paired_with is None
        else:
            assert self.paired_end is not None and self.paired_with is not None

    @classmethod
    def create_from_ENCODE_dict(
        cls,
        file_dict: dict,
        run_type: RunType,
        experiment: str,
        biosample: str,
        technical_replicate_number: int,
        biological_replicate_number: int,
        control: bool,
        antibody: str,
    ):
        "Create ENCODE File object from a dictionary(part of expr_dict)"
        return SE_File(
            file_dict.get("accession"),
            file_dict.get("read_count"),
            file_dict.get("file_format"),
            file_dict.get("no_file_available"),
            file_dict.get("platform", {}).get("term_name"),
            file_dict.get("read_length"),
            experiment,
            control,
            file_dict.get("replicate", {}).get("library")[11:-1],
            biosample,
            technical_replicate_number,
            biological_replicate_number,
            antibody,
            file_dict.get("href"),
            run_type,
            file_dict.get("paired_end", None),
            file_dict.get("paired_with", None),
        )

    def download(self, download_dr: Path | str) -> SE_Fastq:
        """Download the file from ENCODE server to the specified directory. By default downloads zipped files,
        which are saved as .fastq.gz files."""

        if self.no_file_available:
            raise NoFileAvailable(f"No file available for {self.accession}")

        # Constructing the filnames and filepaths
        filename = f"{self.accession}.{self.file_format}"
        filepath = Path(download_dr) / Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # Dont' download if already downloaded
        if filepath.exists() or Path(str(filepath) + ".gz").exists():

            if Path(str(filepath) + ".gz").exists():
                filepath = Path(str(filepath) + ".gz")

            return SE_Fastq(
                filepath,
                self.accession,
                self.platform,
                self.read_length,
                self.experiment,
                self.library,
                self.biosample,
                self.technical_replicate_number,
                self.biological_replicate_number,
                self.control,
                self.antibody,
                self.href,
                self.run_type,
                self.paired_end,
                self.paired_with,
            )

        # Sending request to ENCODE for download stream
        response = requests.get(
            "https://www.encodeproject.org" + self.href, stream=True
        )

        # Check status code of request
        if response.status_code != 200:
            raise FileDownloadError(f"File download error for {self.accession}.")

        # Download the file in chunks
        with filepath.open(mode="wb") as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                file.write(chunk)

        return SE_Fastq(
            filepath,
            self.accession,
            self.platform,
            self.read_length,
            self.experiment,
            self.library,
            self.biosample,
            self.technical_replicate_number,
            self.biological_replicate_number,
            self.control,
            self.antibody,
            self.href,
            self.run_type,
            self.paired_end,
            self.paired_with,
        )

    def update_database(self):
        "Adds information to the Selex_X_Genome.db database"
        # Connect to the database
        with sqlite3.connect(
            "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
        ) as conn:
            cursor = conn.cursor()

            # Get the experiment_id from the experiments table
            try:
                cursor.execute(
                    """
                    SELECT "id" FROM "experiments" WHERE "accession" = ?
                    """,
                    (self.experiment,),
                )
                experiment_id = cursor.fetchone()[0]
            except:
                raise ValueError("Could not find experiment in Database")

            # Get the library_id from the experiments table
            try:
                cursor.execute(
                    """
                    SELECT "id" FROM "libraries" WHERE "accession" = ?
                    """,
                    (self.library,),
                )
                library_id = cursor.fetchone()[0]
            except:
                raise ValueError("Could not find library in Database")

            # Value to be added to DB
            file = [(self.accession, self.read_length, experiment_id, library_id, "n")]
            # Update DB
            try:
                cursor.executemany(
                    """
                INSERT INTO "files" ("accession", "read_length", "experiment_id", "library_id", "paired_end") VALUES (?, ?, ?, ?, ?)
                """,
                    file,
                )
                conn.commit()
            except sqlite3.IntegrityError:
                pass
            except Exception as e:
                print("Error:", e)


class PE_File(ENCODE_Object):
    """Class for a ENCODE file (ENCFFXXXXXX) object"""

    def __init__(self, r1: SE_File, r2: SE_File):
        self.r1 = r1
        self.r2 = r2

    # R1 and R2 are commutative
    def download(self, download_dr: Path) -> PE_Fastq:
        """Download the pe files from ENCODE server to the specified directory."""
        if self.r1.no_file_available or self.r2.no_file_available:
            raise NoFileAvailable(
                f"No file available for {self.r1.accession} and {self.r2.accession}"
            )

        r1_fastq = self.r1.download(download_dr)
        r2_fastq = self.r2.download(download_dr)

        return PE_Fastq(r1_fastq, r2_fastq)
