from pathlib import Path

import pandas as pd
from utils.slurmjob import Slurmjob

df = pd.read_csv(
    "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/cluster_1.csv", index_col=0
)


job_script = Path(
    "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/source/utils/pyprob_2.py"
)

# Slurmjob("/burg/home/hg2604/hblab/Projects/Selex-X-Genome/temp")
for index, row in df.iterrows():

    exp = row["Experiment_Acc.x"]
    lib = row["Library_Acc.x"]
    file = row["File_Acc"] + ".tsv.gz"

    cnt_table = f"/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/CTCF_Human/{exp}/{lib}/{file}"

    slurm_job = Slurmjob(
        Path(f"/burg/home/hg2604/hblab/Projects/Selex-X-Genome/temp/{file}.job"),
        file,
        job_script,
        (cnt_table,),
        file,
        cores=6,
        time=6,
    )
    slurm_job.submitJob()
