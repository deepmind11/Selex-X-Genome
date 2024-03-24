from pathlib import Path

import requests

from .base import ENCODE_Object
from .disk_files import fastq


class File(ENCODE_Object):
    """Class for a ENCODE file (ENCFFXXXXXX) object"""

    def __init__(
        self,
        accession,
        read_count,
        file_format,
        no_file_available,
        platform,
        read_length,
        library,
        href,
        run_type,
        paired_with=None,
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
        self.paired_with: File = None

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
            return fastq(filepath)

        # Downloading the file in chunks
        response = requests.get(self.href, stream=True)
        with filepath.open(mode="wb") as file:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                file.write(chunk)

        return fastq(filepath)


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
