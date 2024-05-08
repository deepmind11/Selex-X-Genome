import subprocess
from pathlib import Path


class DiskFile:
    """Base class for files on disk"""

    def __init__(self, file_path: Path):

        # Path on disk
        self.file_path = file_path
        # Set deleted to False upon initialization
        self.deleted = False
        # Check that file exists upon initialization
        if not self.file_path.exists():
            raise FileNotFoundError(f"File {self.file_path} does not exist.")
        # Set flag for zipped
        self.zipped = DiskFile.is_gz_file(self.file_path)
        # If file is zipped make sure it has a .gz extension
        if self.zipped and self.file_path.suffix != ".gz":
            self.file_path.rename(str(self.file_path) + ".gz")
            self.file_path = Path(str(self.file_path) + ".gz")

    @staticmethod
    def is_gz_file(filepath):
        with open(filepath, "rb") as test_f:
            return test_f.read(2) == b"\x1f\x8b"

    def delete(self):
        if self.file_path is not None:
            self.file_path.unlink()
            self.file_path = None
            self.deleted = True

    def zip(self):
        "Zip file"
        if not self.deleted and not self.zipped:
            subprocess.run(["gzip", str(self.file_path)])
            # Update path
            self.file_path = Path(str(self.file_path) + ".gz")
            # Update zipped flag
            self.zipped = DiskFile.is_gz_file(self.file_path)

    def unzip(self):
        "Unzip File"
        if not self.deleted and self.zipped:
            subprocess.run(["gzip", "-d", str(self.file_path)])
            # Update path
            self.file_path = Path(str(self.file_path)[:-3])
            # Update zipped flag
            self.zipped = DiskFile.is_gz_file(self.file_path)
