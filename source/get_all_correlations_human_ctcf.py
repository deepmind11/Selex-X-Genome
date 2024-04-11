import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.append("/burg/hblab/users/hg2604/Projects/OOP/Selex-X-Genome/source")
from source.diskfiles.slurmjob import Slurmjob

list_of_correlations = list(
    Path("/burg/home/hg2604/hblab/Projects/OOP/Selex-X-Genome/data/CTCF_Human").rglob(
        "*correlation.txt"
    )
)


dict_list = []

for file in list_of_correlations:

    row = {}
    row["File_Acc"] = file.name[:11]
    row["Library_Acc"] = file.parent.name
    row["Experiment_Acc"] = file.parent.parent.name

    with open(file, "r") as f:
        row["Correlation"] = float(f.readline().strip("\n"))

    dict_list.append(row)


CORRELATION_DF = pd.DataFrame(dict_list)

CORRELATION_DF.to_csv(
    Path(
        "/burg/home/hg2604/hblab/Projects/OOP/Selex-X-Genome/CTCF_MOTIF_CORRELATION.csv"
    ),
    index=False,
)
