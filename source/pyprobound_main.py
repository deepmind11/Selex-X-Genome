"""Create the scatter plot and correlation for all the count tables"""

import subprocess
import sys
from pathlib import Path

sys.path.append("/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source")
from source.diskfiles.slurmjob import Slurmjob

list_of_tables = list(
    Path("/burg/home/hg2604/hblab/Projects/OOP/Selex-X-Genome/data/CTCF_Human").rglob(
        "*.tsv.gz"
    )
)


for table in list_of_tables:
    # print(table)
    # Create a Slurmjob object
    jobfile_Path = table.parent / Path(table.stem + ".job")
    job_name = table.stem
    job_script = Path(
        "/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source/pyprobound_script.py"
    )
    job_params = (table,)

    sj = Slurmjob(jobfile_Path, job_name, job_script, job_params)
    sj.create_file()
    sj.submitJob()

    # break
