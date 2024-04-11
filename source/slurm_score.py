#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path

sys.path.append("/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source")

from base import DiskFile, ENCODE_Object, Experiment
from experiment import Control, TFChipSeq
from search import EncodeSearch

from source.diskfiles.slurmjob import Slurmjob
from source.motifs.motif import Mononucleotide
from source.motifs.parse_motifcentral_json import MOTIFCENTRAL

max_motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[222])


Human_MYC = EncodeSearch("CTCF", "Homo+sapiens", limit="all")
Human_MYC.search()

# print(Human_MYC.search_result)
Experiments = Human_MYC.get_experiments()

for experiment in Experiments:

    slurmjob_path = Path(
        f"/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/data/CTCF_Human/{experiment.accession}/score.job"
    )
    job_name = "Score_count_tables"
    job_script = Path(
        "/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source/score_experiment_cnt_tbl.py"
    )
    job_params = (experiment.accession,)

    # Create slurm job file
    job_file = Slurmjob(
        file_path=slurmjob_path,
        job_name=job_name,
        job_script=job_script,
        job_params=job_params,
    )
    job_file.create_file()
    job_file.submitJob()
