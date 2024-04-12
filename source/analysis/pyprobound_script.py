#!/usr/bin/env python3

import os
import sys
from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyprobound
import pyprobound.plotting
import torch

parser = ArgumentParser()
parser.add_argument("cntTblPath", type=str)


args = parser.parse_args()


# ! Take count table as input
count_table_path = Path(args.cntTblPath)


alphabet = pyprobound.alphabets.DNA()


# In[3]:
dataframe = pd.read_csv(
    count_table_path,
    header=None,
    index_col=0,
    sep="\t",
)


# In[4]:
count_table = pyprobound.CountTable(dataframe, alphabet, zero_pad=True)


# In[6]:
nonspecific = pyprobound.layers.NonSpecific(alphabet=alphabet, name="NS")
psams = [
    pyprobound.layers.PSAM(
        alphabet=alphabet,
        kernel_size=18,
        seed=["---CGCCMYCTAGTGG--"],
        name="CTCF",
    )
]


# In[7]:

conv0d = pyprobound.layers.Conv0d.from_nonspecific(nonspecific, count_table)
conv1ds = [
    pyprobound.layers.Conv1d.from_psam(
        psam, count_table, train_posbias=True, length_specific_bias=False, bias_bin=5
    )
    for psam in psams
]

modes = [pyprobound.Mode([conv0d])] + [pyprobound.Mode([conv1d]) for conv1d in conv1ds]

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
    checkpoint=count_table_path.parent / Path(count_table_path.stem + ".pt"),
    output=count_table_path.parent / Path(count_table_path.stem + ".txt"),
)


optimizer.train_sequential()
optimizer.reload()

with torch.inference_mode():
    loss, reg = model([count_table])
    print(loss, reg, loss + reg)

# Generating all the plots

# Generating all the plots
for index, psam in enumerate(psams):
    logo_for = pyprobound.plotting.logo(psam)
    logo_for.savefig(
        count_table_path.parent / Path(count_table_path.stem + "_forwardLogo.png")
    )
    logo_rev = pyprobound.plotting.logo(psam, reverse=True)
    logo_rev.savefig(
        count_table_path.parent / Path(count_table_path.stem + "_reverseLogo.png")
    )

# Pos Bias Profile
for index, conv1d in enumerate(conv1ds):
    pos_bias = pyprobound.plotting.posbias(conv1d)
    pos_bias.savefig(
        count_table_path.parent / Path(count_table_path.stem + "_positionBias.png")
    )

# #Probe Enrichment
probe_er = pyprobound.plotting.probe_enrichment(experiment, count_table)
probe_er.savefig(
    count_table_path.parent / Path(count_table_path.stem + "_CTCF_probeEnrichment.png")
)

# #Mode Contribution
mode_cr = pyprobound.plotting.contribution(round_1, count_table)
mode_cr.savefig(
    count_table_path.parent / Path(count_table_path.stem + "_CTCF_modeContribution.png")
)


# Compute the correlation between
CTCF_SELEX = [
    [
        -0.3238,
        0.2350,
        -0.3666,
        -0.7458,
        -0.0185,
        -0.4299,
        -0.5675,
        1.2223,
        -0.7196,
        0.1198,
        -1.2056,
        -0.9507,
        -1.6197,
        -0.4313,
        -0.1272,
        -0.6833,
        -0.3878,
        0.0664,
    ],
    [
        -0.2858,
        -0.8964,
        0.9162,
        1.5294,
        -0.8656,
        0.6523,
        -0.2411,
        -1.1736,
        -1.2920,
        -1.4794,
        -0.8737,
        -1.0614,
        -0.6112,
        1.3176,
        -0.7485,
        0.5472,
        -0.2614,
        -0.1690,
    ],
    [
        0.0152,
        0.0253,
        -0.8351,
        -0.5783,
        0.1412,
        -0.2020,
        -0.8803,
        -0.5201,
        1.8074,
        1.1590,
        0.9481,
        1.8797,
        1.4305,
        -0.9774,
        1.0567,
        -0.6921,
        -0.8604,
        -0.2929,
    ],
    [
        -0.2603,
        -0.6660,
        -0.8142,
        -1.3049,
        -0.3568,
        -1.1201,
        0.5892,
        -0.6283,
        -0.8955,
        -0.8991,
        0.0316,
        -0.9674,
        -0.2993,
        -1.0086,
        -1.2806,
        -0.3897,
        0.3424,
        -0.6921,
    ],
]


current_psam = psams[0].get_filter(0)[1, :, :].tolist()


def compare_psam_scatterplot(
    psam1, psam2, xlabel="SELEX-CTCF", ylabel="ESTIMATED-CTCF"
):
    # The psams need to be the same length
    # Generate random data for the scatter plots

    x1 = psam1[0]
    x2 = psam1[1]
    x3 = psam1[2]
    x4 = psam1[3]
    y1 = psam2[0]
    y2 = psam2[1]
    y3 = psam2[2]
    y4 = psam2[3]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the scatter plots with different colors
    ax.scatter(x1, y1, color="green", label="A")
    ax.scatter(x2, y2, color="blue", label="C")
    ax.scatter(x3, y3, color="yellow", label="G")
    ax.scatter(x4, y4, color="red", label="T")

    # Add a legend for the colors
    ax.legend()

    # Set labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("PSAM Scatter Plot")

    # Display the plot
    plt.savefig(
        count_table_path.parent / Path(count_table_path.stem + "_PSAM_scatter_plot.png")
    )


compare_psam_scatterplot(CTCF_SELEX, current_psam)


def get_correlation(psam1, psam2):
    """Returns the pearson correalation between psam1 and psam2"""

    # The psams need to be the same length
    psam1 = psam1[0] + psam1[1] + psam1[2] + psam1[3]
    psam2 = psam2[0] + psam2[1] + psam2[2] + psam2[3]

    psam1 = np.array(psam1, dtype="float32")
    psam2 = np.array(psam2, dtype="float32")

    correlation = np.corrcoef(psam1, psam2)
    correlation = correlation[0, 1]

    # Saving the correlation to file:
    correlation_path = count_table_path.parent / Path(
        count_table_path.stem + "_correlation.txt"
    )
    correlation_path.touch()
    with open(correlation_path, "w") as f:
        f.write(str(correlation))

    return


get_correlation(CTCF_SELEX, current_psam)
