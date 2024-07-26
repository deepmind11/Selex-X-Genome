import math
from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":

    parser = ArgumentParser()

    # Required positional arguments
    parser.add_argument("CountTablePath", help="Absolute Path to the Count Table")
    parser.add_argument("Motif_ID", help="Motif_ID")

    args = parser.parse_args()

    CountTablePath = Path(args.CountTablePath)
    Motif_ID = args.Motif_ID

    if not CountTablePath.exists():
        raise FileNotFoundError("The counttable does not exist")

    ScorePath = CountTablePath.parent / Path(
        f"{CountTablePath.name[:11]}_{Motif_ID}.txt"
    )

    if not ScorePath.exists():
        raise FileNotFoundError(
            f"You need to score the count table with the motif first!"
        )

    # Reading data as pandas dataframe
    pandas_r0_r1 = pd.read_csv(
        CountTablePath, sep="\t", header=None, index_col=None, usecols=[1, 2]
    )
    score_pd = pd.read_csv(ScorePath, sep="\t", header=None, index_col=None)

    merged_df = pd.concat([pandas_r0_r1, score_pd], axis=1)
    merged_df.columns = ["R0", "R1", "Score"]
    merged_df.sort_values(by="Score", ascending=False, inplace=True)

    merged_df = merged_df.apply(pd.to_numeric)
    merged_df.reset_index(drop=True, inplace=True)

    # Enrirchment vs Bin Size plots
    num_rows = merged_df.shape[0]
    bin_sizes = range(500, num_rows)
    enrichment = []

    numerator = merged_df.loc[: 500 - 1, "R1"].sum() + 1
    denominator = merged_df.loc[: 500 - 1, "R0"].sum() + 1

    enrichment.append(numerator / denominator)

    for index, row in merged_df.loc[501:].iterrows():
        numerator += row["R1"]
        denominator += row["R0"]
        enrichment.append(numerator / denominator)

    # Creating the scatter plot
    plt.scatter([math.log(bin, 10) for bin in bin_sizes], enrichment)

    # Adding titles and labels
    plt.title(f"Enrichment vs Bin-Size for {CountTablePath.name[:11]} using {Motif_ID}")
    plt.xlabel("Log(10) Bin-Size")
    plt.ylabel("Enrichment for Top Bin")

    # Saving the plot to a specific location
    save_path = CountTablePath.parent / Path(
        f"{CountTablePath.name[:11]}_{Motif_ID}_Enrichment_vs_BinSize.png"
    )
    plt.savefig(save_path)
