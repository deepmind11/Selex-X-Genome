from pathlib import Path

import matplotlib.pyplot as plt
from diskfiles.countTables import CountTable

table = CountTable(
    Path(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/data/CTCF_Homo_sapiens/ENCSR000AUE/ENCLB695AKO/ENCFF000AHR.tsv.gz"
    )
)

from motifs.parse_motifcentral_json import MOTIFCENTRAL

gene = "CTCF"
for i, motif in enumerate(MOTIFCENTRAL):
    if gene in motif["metadata"]["factors"][0]["gene_symbol"]:
        print(i)

from motifs.motif import Mononucleotide

motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[222])

plot = table.plot_enrichment_vs_bin(motif)

fig = plot.get_figure()
fig.savefig("/burg/home/hg2604/hblab/Projects/Selex-X-Genome/plot.png")
