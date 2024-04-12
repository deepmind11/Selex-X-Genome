#!/usr/bin/env python3

import multiprocessing
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyprobound
import pyprobound.plotting
import torch

dr_path = Path("/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/CTCF_Human")
pt_file = list(dr_path.rglob("*.tsv.pt"))
txt_file = list(dr_path.rglob("*.tsv.txt"))

input_tuples = list(zip(pt_file, txt_file))


# Calling the class as a function...this is enabled by something from the typing module
alphabet = pyprobound.alphabets.DNA()

# In[3]:
dataframe = pd.read_csv(
    "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/CTCF_Human/ENCSR000AHD/ENCLB305ZZZ/ENCFF000QLW.tsv.gz",
    header=None,
    index_col=0,
    sep="\t",
)

# In[4]:
count_table = pyprobound.CountTable(dataframe, alphabet, zero_pad=True)


def add_psam_to_df(input_tuple):

    pt = str(input_tuple[0])
    txt_file = str(input_tuple[1])

    # In[6]:
    nonspecific = pyprobound.layers.NonSpecific(alphabet=alphabet, name="NS")
    psams = [
        pyprobound.layers.PSAM(
            alphabet=alphabet,
            kernel_size=18,
            seed=["---CGCCMYCTAGTGG--"],
            name="CTCF",
            shift_footprint_heuristic=True,
            increment_footprint=True,
        )
    ]

    conv0d = pyprobound.layers.Conv0d.from_nonspecific(nonspecific, count_table)
    conv1ds = [
        pyprobound.layers.Conv1d.from_psam(
            psam,
            count_table,
            train_posbias=True,
            length_specific_bias=False,
            bias_bin=5,
        )
        for psam in psams
    ]

    modes = [pyprobound.Mode([conv0d])] + [
        pyprobound.Mode([conv1d]) for conv1d in conv1ds
    ]

    round_0 = pyprobound.rounds.InitialRound()
    round_1 = pyprobound.rounds.BoundUnsaturatedRound.from_binding(modes, round_0)

    experiment = pyprobound.Experiment(
        [round_0, round_1],
        name="CTCF",
        counts_per_round=count_table.counts_per_round,
    )

    model = pyprobound.MultiExperimentLoss([experiment], pseudocount=20)

    optimizer = pyprobound.Optimizer(
        model,
        [count_table],
        greedy_threshold=2e-4,
        device="cpu",
        checkpoint=pt,
        output=txt_file,
    )

    optimizer.reload()

    del_G = (
        list(psams[0].get_filter(0)[1, :, :].cpu().detach().numpy()[0])
        + list(psams[0].get_filter(0)[1, :, :].cpu().detach().numpy()[1])
        + list(psams[0].get_filter(0)[1, :, :].cpu().detach().numpy()[2])
    )
    lib = Path(pt).parent.stem
    exp = Path(pt).parent.parent.stem
    file = Path(pt).name[:11]

    return [exp, lib, file] + del_G


with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
    # Map the function to the list of items and collect the results
    results = pool.map(add_psam_to_df, input_tuples)

filt_results = list(filter(lambda x: len(x) == 57, results))


colnames_x = ["Exp", "Lib", "File"]
for j in ["A", "C", "G"]:
    for i in range(18):
        colnames_x.append(j + str(i))


df_x = pd.DataFrame(filt_results, columns=colnames_x)

df_x.to_csv(
    "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/one_off_scripts/psam_df.csv",
    index=False,
)
