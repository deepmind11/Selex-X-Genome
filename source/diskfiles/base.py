import subprocess
from pathlib import Path


class DiskFile:
    """Base class for files on disk"""

    def __init__(
        self, file_path: Path = None
    ):  # Polymorphism in the tpyes supported by file_path. For ex: Path | list(Path)
        self.file_path = file_path

    def delete(self):
        if self.file_path is not None:
            self.file_path.unlink()
            self.file_path = None

    def zip(self):
        if self.file_path is not None:
            if self.file_path.suffix != ".gz":
                subprocess.run(["gzip", str(self.file_path)])
                # Update path
                self.file_path = Path(str(self.file_path) + ".gz")
