import subprocess
import sys
from pathlib import Path

from diskfiles.base import DiskFile


class Slurmjob(DiskFile):
    """Class for creating a Slurmjob files."""

    def __init__(
        self,
        file_path: Path,
        job_name: str,
        job_script: Path,
        job_params: tuple,
        output: str,
        cores: int = 4,
        time: int = 4,
    ):
        self.file_path = file_path
        self.job_name = job_name
        self.job_script = job_script  # has to be a script py, sh, anything
        self.job_params = job_params
        self.output = output
        self.cores = cores
        self.time = time
        self.create_file()
        super().__init__(file_path)

    def create_file(self):
        self.file_path.parent.mkdir(exist_ok=True, parents=True)
        self.file_path.touch(exist_ok=True)
        with self.file_path.open(mode="w") as jf:
            jf.write("#!/bin/bash\n")
            jf.writelines(f"#SBATCH --job-name={self.job_name}\n")
            jf.writelines(
                f"#SBATCH --output={self.file_path.parent}/{self.output}.out\n"
            )
            jf.writelines(
                f"#SBATCH --error={self.file_path.parent}/{self.output}.err\n"
            )
            jf.writelines(f"#SBATCH -c {self.cores}\n")
            jf.writelines("#SBATCH --mem-per-cpu=5G\n")
            jf.writelines("#SBATCH --account=hblab\n")
            jf.writelines(f"#SBATCH -t {self.time}:59:00\n\n")
            jf.writelines(
                f'{str(self.job_script)} {" ".join(str(param) for param in self.job_params)}'
            )

    def submitJob(self):
        if not self.file_path.exists():
            self.create_file()

        # Call subprocess
        subprocess.run(["sbatch", str(self.file_path)])
