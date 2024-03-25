from pathlib import Path

import requests

from .base import ENCODE_Object
from .disk_files import PE_Fastq, SE_Fastq


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
        paired_with=None,  # WILL be none for SE (Assert this)
    ):
        self.accession = accession
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

        # file_to_download["paired_with"][7:-1]

        # How do positional and optional arguments work in python
        # A constructor method that creates a encode frile from the dict returned by REST API

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
        response = requests.get(self.href, stream=True)
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

    # How do positional and optional arguments work in python
    # A constructor method that creates a encode frile from the dict returned by REST API

    # Order does not matter as long as the pairs are in the same folder
    def download(self, download_dr: Path):
        """Download the pe files from ENCODE server to the specified directory."""
        r1_fastq = self.r1.download(download_dr)
        r2_fastq = self.r2.download(download_dr)

        return PE_Fastq(r1_fastq, r2_fastq)


#     'accession', str
# 2. 'read_count', int **Verify its over 2 million**
# 3. 'file_format', str **Super Important**
# 4. 'no_file_available', bool
# 5. 'platform', dict
# 6. 'read_length', int

# 8.
# 9. ['replicate']['library'] **Super Important**
# 10. 'href', str   **This is the downlaod link**

# 13. ['run_type'] **Important to distinguish b
