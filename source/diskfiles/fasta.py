from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from diskfiles.countTables import CountTable
    from ENCODE.base import RunType

from diskfiles.base import DiskFile
from ENCODE.experiment import TFChipSeq


class CantBuildTableFromControl(Exception):
    pass


class AmbiguousControls(Exception):
    pass


class SE_Fasta(DiskFile):
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

    def build_count_table(self) -> CountTable:
        """Finds the matching control and builds the count table.(Returns Zipped)"""
        # Establish the R0 and R1 Fasta files
        if self.control:
            raise CantBuildTableFromControl("Call from ChIPseq experiment instead.")
        else:
            # R1
            r1 = self.file_path
            # Fetching the appropriate controls
            exp = TFChipSeq(self.experiment)
            control = exp.get_controls()
            control_libs = control.get_libraries()
            # Filtering the libraries based on Biosample
            control_libs = list(
                filter(
                    lambda x: (x.biosample == self.biosample),
                    control_libs,
                )
            )
            # If more than one control lib then do further filtering based on technical replicate number
            if len(control_libs) > 1:
                control_libs = list(
                    filter(
                        lambda x: (
                            x.technical_replicate_number
                            == self.technical_replicate_number
                        ),
                        control_libs,
                    )
                )
            # If length controls libs not equal 1, raise error
            if len(control_libs) != 1:
                raise AmbiguousControls
            # Get the files for the control lib
            control_file = control_libs[0].get_Files()[0]
            # Path to the hypothetical fasta file
            r0 = Path(
                f"/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/Control/{control.accession}/{control_libs[0].accession}/{control_file.accession}.fasta"
            )
            # Creating the count table.
            cnt_tbl_path = self.file_path.parent / Path(f"{self.accession}.tsv")
            if r0.exists():
                count_table = CountTable.create_from_fasta(
                    r0, r1, cnt_tbl_path
                )  # Returns if already exists
                count_table.zip()
                return count_table
            else:
                # Downlaoding all the required files
                control_fastq = control_file.download(r0.parent)
                control_fastq = (
                    control_fastq.subsample()
                )  #! Currently happening with the defaults.
                control_fasta = control_fastq.transform_to_fasta()
                assert r0 == control_fasta.file_path
                # Returning the count table
                count_table = CountTable.create_from_fasta(r0, r1, cnt_tbl_path)
                count_table.zip()
                return count_table
