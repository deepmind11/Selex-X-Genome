from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":

    parser = ArgumentParser()

    # Required positional arguments
    parser.add_argument("CountTablePath", help="Absolute Path to the Count Table")
    parser.add_argument("Ref_Motif_ID", help="Ref_Motif_ID")
    parser.add_argument("Cofactor_Motif_ID", help="Cofactor_Motif_ID")

    args = parser.parse_args()

    CountTablePath = Path(args.CountTablePath)
    Ref_Motif_ID = args.Ref_Motif_ID
    Cofactor_Motif_ID = args.Cofactor_Motif_ID

    if not CountTablePath.exists():
        raise FileNotFoundError("The counttable does not exist")

    Ref_ScorePath = CountTablePath.parent / Path(
        f"{CountTablePath.name[:11]}_{Ref_Motif_ID}.txt"
    )

    if not Ref_ScorePath.exists():
        raise FileNotFoundError(
            f"You need to score the count table with the reference motif first!"
        )

    Cofactor_ScorePath = CountTablePath.parent / Path(
        f"{CountTablePath.name[:11]}_{Cofactor_Motif_ID}.txt"
    )

    if not Cofactor_ScorePath.exists():
        raise FileNotFoundError(
            f"You need to score the count table with the cofactor motif first!"
        )

    # Reading data as pandas dataframe
    pandas_r0_r1 = pd.read_csv(
        CountTablePath, sep="\t", header=None, index_col=None, usecols=[1, 2]
    )
    ref_score_pd = pd.read_csv(Ref_ScorePath, sep="\t", header=None, index_col=None)
    cofactor_score_pd = pd.read_csv(
        Cofactor_ScorePath, sep="\t", header=None, index_col=None
    )

    merged_df = pd.concat([pandas_r0_r1, ref_score_pd, cofactor_score_pd], axis=1)
    merged_df.columns = ["R0", "R1", "Ref_Score", "Cofactor Score"]
    merged_df.sort_values(by="Ref_Score", ascending=False, inplace=True)

    # Bin the dataframe

    i = 0
    avg_ref_score = []
    avg_cof_score = []
    r0_count = []
    r1_count = []
    while (
        i < merged_df.shape[0] - 10
    ):  # Keeping a buffer, so that the last bin isn't super small
        bin = merged_df.iloc[i : i + 1000, :].copy()
        avg_cof_score.append(bin.iloc[:, 3].mean())
        avg_ref_score.append(bin.iloc[:, 2].mean())
        r0_count.append(bin.iloc[:, 0].sum())
        r1_count.append(bin.iloc[:, 1].sum())
        i += 1000

    bin_df = pd.DataFrame(
        {
            "avg_cof_score": avg_cof_score,
            "avg_ref_score": avg_ref_score,
            "r0_count": r0_count,
            "r1_count": r1_count,
        }
    )
    bin_df["enrichment"] = bin_df["r1_count"] / bin_df["r0_count"]

    # Enrirchment vs TF and Co-factor Scores
    num_rows = list(range(bin_df.shape[0]))

    enrichment = bin_df["enrichment"]
    avg_cof_score = bin_df["avg_cof_score"]
    avg_ref_score = bin_df["avg_ref_score"]

    # Creating the scatter plot
    # Plotting the data
    plt.plot(num_rows, enrichment, label="Enrichment", color="r")  # red color for y1
    plt.plot(
        num_rows, avg_cof_score, label=Cofactor_Motif_ID, color="g"
    )  # green color for y2
    plt.plot(
        num_rows, avg_ref_score, label=Ref_Motif_ID, color="b"
    )  # blue color for y3

    # Adding a legend
    plt.legend()

    # Adding titles and labels
    plt.title("Enrichment, Score, Co-factor score vs Bin Number(sorted by Score)")
    plt.xlabel("Bin Number")
    plt.ylabel("Y-axis")

    # Saving the plot to a specific location
    save_path = CountTablePath.parent / Path(
        f"{CountTablePath.name[:11]}_{Ref_Motif_ID}_{Cofactor_Motif_ID}.png"
    )
    plt.savefig(save_path)
