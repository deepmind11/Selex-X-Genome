from pathlib import Path


class Slurmjob:
    """Class for creating a Slurmjob files."""

    def __init__(self, job_name, job_description, job_script, job_output):
        self.job_name = job_name
        self.job_description = job_description
        self.job_script = job_script
        self.job_output = job_output

    def create(self) -> Path:
        pass
