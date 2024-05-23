from __future__ import annotations

import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

from diskfiles.base import DiskFile
from ENCODE.base import RunType
from utils.slurmjob import Slurmjob

if TYPE_CHECKING:
    from diskfiles.fasta import SE_Fasta

# File size (number of reads): implement methods

# Other files in the library(series)
# Possible control fastq files (ideally just one) (this will be used to make the count table)
# The ability to make count_tables should be a method from here. Instead of orchestrating everything externally.

# Need to clearly establish the file structure for all my downloads. (Search for library folder in the data folder)


class SE_Fastq(DiskFile):
    """Class for SE fastq files on disk"""

    def __init__(
        self,
        file_path: Path,
        accession: str,
        platform: str,
        read_length: int,
        experiment: str,
        library: str,
        biosample: str,
        technical_replicate_number: int,
        biological_replicate_number: int,
        control: bool,
        antibody: str,
        href: str,
        run_type: RunType,  # This is either SE or PE
        paired_end=None,  # Will be none for SE (Assert this)
        paired_with: str = None,
    ):
        super().__init__(file_path)
        self.accession = accession
        self.platform = platform
        self.read_length = read_length
        self.experiment = experiment
        self.library = library
        self.biosample = biosample
        self.technical_replicate_number = technical_replicate_number
        self.biological_replicate_number = biological_replicate_number
        self.control = control
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

    def subsample(self, seed: int = 42, size: int = 1000000) -> SE_Fastq:
        """Returns a new subsampled SE_Fastq object"""
        # Create path for subsampled fastq
        subsampled_fastq_file = self.file_path.parent / Path(
            self.accession + "_subsampled.fq"
        )

        # If subsampled_fastq_file exists, then return it
        if subsampled_fastq_file.exists():
            return SE_Fastq(
                subsampled_fastq_file,
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

        # Otherwise, call seqtk to subsample
        else:
            with open(subsampled_fastq_file, "w") as output_file:
                subprocess.run(
                    ["seqtk", "sample", f"-s{seed}", str(self.file_path), str(size)],
                    stdout=output_file,
                )

            return SE_Fastq(
                subsampled_fastq_file,
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

    def transform_to_fasta(
        self,
        transform_script: Path = Path(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/source/se_processing_template.sh"
        ),
    ) -> SE_Fasta:
        """Tranfrom fastq to fasta for now. FileAcc.fasta"""
        from diskfiles.fasta import SE_Fasta

        # If exists then return
        fasta_file = self.file_path.parent / Path(self.accession + ".fasta")
        if fasta_file.exists():
            return SE_Fasta(
                fasta_file,
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
        else:
            # In my use case, it should produce a fasta file in the same directory.
            subprocess.run([str(transform_script), str(self.file_path)])
            return SE_Fasta(
                fasta_file,
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

    def slurm_transform_to_fasta(  # ! Need to implement this properly. Wait for the slurm job to finish.
        self,
        transform_script: Path = Path(
            "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/source/se_processing_template.sh"
        ),
    ) -> SE_Fasta:
        # Create slurm job file
        fasta_file = self.file_path.parent / Path(self.accession + ".fasta")
        if fasta_file.exists():
            return
        else:
            slurm_file_path = self.file_path.parent / Path(
                f"{self.accession}_transfrom.job"
            )
            slurm_job = Slurmjob(
                file_path=slurm_file_path,
                job_name=f"{self.accession}_transfrom",
                job_script=transform_script,
                job_params=(str(self.file_path),),
                output=f"{self.accession}_transfrom",
            )
            slurm_job.create_file()
            slurm_job.submitJob()

    def download(self):
        """Downloads file again if deleted."""
        pass

    def get_all_files_in_library(self) -> list["SE_Fastq"]:
        """Returns a list of all fastq files belonging to this library."""
        pass

    def get_control_files(self) -> list["SE_Fastq"]:
        """Returns a list of all control files for this library."""
        pass


# This can be implement as a tuple of SE_Fastq objects.
class PE_Fastq:
    """Class for PE fastq files on disk"""

    def __init__(self, r1: SE_Fastq, r2: SE_Fastq):
        self.r1 = r1
        self.r2 = r2

        # ! Assert that both are None or not None together.
        assert self.r1.run_type == RunType.PE and self.r2.run_type == RunType.PE

    def delete(self) -> None:
        """Deletes both the fastq files from disk."""
        self.r1.delete()
        self.r2.delete()

    def subsample(self, seed=42, size=1000000):  # For now let's keep seed 42
        """Returns subsampled fastq object"""

        subsampled_read1 = self.r1.subsample(seed=seed, size=size)
        subsampled_read2 = self.r2.subsample(seed=seed, size=size)

        return PE_Fastq(subsampled_read1, subsampled_read2)

    def transform(self, transform_script: Path):
        # In my use case, it should produce a fasta file in the same directory.
        subprocess.run([transform_script, self.r1.file_path, self.r2.file_path])

    def download(self):
        """Downloads file again if deleted."""
        pass
