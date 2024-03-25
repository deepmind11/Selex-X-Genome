import subprocess
import sys
from pathlib import Path

sys.path.append("/Users/hgz/Documents/researchHBLab/Projects/Selex-X-Genome/source")

from base import DiskFile


class SE_Fastq(DiskFile):
    """Class for SE fastq files on disk"""

    def __init__(self, file_path: Path = None) -> None:
        super().__init__(file_path)

    def subsample(self, seed=42, size=2000000):  # For now let's keep seed 42
        """Returns subsampled fastq object"""

        # Check if file exists
        if not self.file_path.exists():
            raise FileNotFoundError(f"File {self.file_path} does not exist.")

        # Create path for subsampled fastq
        subsampled_fastq_file = self.file_path.parent / Path(
            self.file_path.name[:-6] + "_subsampled.fq"
        )

        # If subsampled_fastq_file exists, then return it
        if subsampled_fastq_file.exists():
            return SE_Fastq(subsampled_fastq_file)

        # Otherwise, call seqtk to subsample
        else:
            with open(subsampled_fastq_file, "w") as output_file:
                subprocess.run(
                    ["seqtk", "sample", f"-s{seed}", str(self.file_path), str(size)],
                    stdout=output_file,
                )

            return SE_Fastq(subsampled_fastq_file)

    def transform(self, transform_script: Path):
        # In my use case, it should produce a fasta file in the same directory.
        subprocess.run([transform_script, self.file_path])


# This can be implement as a tuple of SE_Fastq objects.
class PE_Fastq:
    """Class for PE fastq files on disk"""

    def __init__(self, r1: SE_Fastq = None, r2: SE_Fastq = None) -> None:
        self.r1 = r1
        self.r2 = r2

        # ! Assert that both are None or not None together.

    def delete(self) -> None:
        """Deletes both the fastq files from disk."""
        self.r1.delete()
        self.r2.delete()

    def subsample(self, seed=42, size=2000000):  # For now let's keep seed 42
        """Returns subsampled fastq object"""

        subsampled_read1 = self.r1.subsample(seed=seed, size=size)
        subsampled_read2 = self.r2.subsample(seed=seed, size=size)

        return PE_Fastq(subsampled_read1, subsampled_read2)

    def transform(self, transform_script: Path):
        # In my use case, it should produce a fasta file in the same directory.
        subprocess.run([transform_script, self.r1.file_path, self.r2.file_path])
