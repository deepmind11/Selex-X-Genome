from pathlib import Path

import requests
from diskfiles.fastq import PE_Fastq, SE_Fastq
from ENCODE.base import ENCODE_Object


class SE_File(ENCODE_Object):
    """Class for a SE ENCODE file (ENCFFXXXXXX) object"""

    def __init__(
        self,
        accession,
        read_count,
        file_format,
        no_file_available: bool,
        platform,
        read_length,
        library,  # Belongs to
        href,
        run_type,  # This has two optional parameters pe or se
        paired_end=None,  # Will be none for SE (Assert this)
        paired_with: str = None,  # WILL be none for SE (Assert this)
    ):
        super().__init__(accession)
        self.read_count = read_count
        self.file_format = file_format
        self.no_file_available = no_file_available
        self.platform = platform
        self.read_length = read_length
        self.library = library
        self.href = href
        self.run_type = run_type
        self.paired_end = paired_end
        self.paired_with = paired_with

        # Assert
        if self.run_type == "single-ended":
            assert self.paired_end is None and self.paired_with is None
        else:
            assert self.paired_end is not None and self.paired_with is not None

    @classmethod
    def create_from_ENCODE_dict(cls, file_dict: dict):

        return SE_File(
            file_dict.get("accession"),
            file_dict.get("read_count"),
            file_dict.get("file_format"),
            file_dict.get("no_file_available"),
            file_dict.get("platform", {}).get("term_name"),
            file_dict.get("read_length"),
            file_dict.get("replicate", {}).get("library")[11:-1],
            file_dict.get("href"),
            file_dict.get("run_type"),
            file_dict.get("paired_end", None),
            file_dict.get("paired_with", None),
        )

    def download(self, download_dr: Path):
        """Download the file from ENCODE server to the specified directory."""

        # Constructing the filnames and filepaths
        filename = f"{self.accession}.{self.file_format}"
        filepath = download_dr / Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        # Dont' download if already downloaded
        if filepath.exists():
            return SE_Fastq(filepath)

        # Downloading the file in chunks
        response = requests.get(
            "https://www.encodeproject.org" + self.href, stream=True
        )
        with filepath.open(mode="wb") as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                file.write(chunk)

        return SE_Fastq(filepath)


class PE_File(ENCODE_Object):
    """Class for a ENCODE file (ENCFFXXXXXX) object"""

    def __init__(self, r1: SE_File, r2: SE_File):
        self.r1 = r1
        self.r2 = r2

    # file_to_download["paired_with"][7:-1]

    # R1 and R2 are commutative
    def download(self, download_dr: Path):
        """Download the pe files from ENCODE server to the specified directory."""
        r1_fastq = self.r1.download(download_dr)
        r2_fastq = self.r2.download(download_dr)

        return PE_Fastq(r1_fastq, r2_fastq)
