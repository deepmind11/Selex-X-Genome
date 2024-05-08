from pathlib import Path

import pandas as pd

df = pd.read_csv(
    "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/cluster_1.csv", index_col=0
)

corr_motif1 = []
corr_motif2 = []

for index, row in df.iterrows():

    exp = row["Experiment_Acc.x"]
    lib = row["Library_Acc.x"]
    corr_file = row["File_Acc"] + ".tsv_correlation.txt"

    corr_file_path = Path(
        f"/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/CTCF_Human/{exp}/{lib}/{corr_file}"
    )

    with corr_file_path.open(mode="r") as f:
        for line in f.readlines():
            line = line.strip().split(".")

            try:
                if line[1][-2] == "-":
                    corr_motif1.append(-1 * float("." + line[2][:6]))
                else:
                    corr_motif1.append(float("." + line[2][:6]))

                if line[2][-2] == "-":
                    corr_motif2.append(-1 * float("." + line[3][:6]))
                else:
                    corr_motif2.append(float("." + line[3][:6]))
            except:
                corr_motif1.append(0)
                corr_motif2.append(0)


df["corr_motif1"] = corr_motif1
df["corr_motif2"] = corr_motif2

df.to_csv("/burg/home/hg2604/hblab/Projects/Selex-X-Genome/cluster_1_with_corr.csv")
