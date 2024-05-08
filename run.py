from ENCODE.search import EncodeSearch
from motifs.motif import Mononucleotide
from motifs.parse_motifcentral_json import MOTIFCENTRAL

nr3c1 = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[867])

nr3c1.motif_X_TF(EncodeSearch("NR3C1", "Homo+sapiens"))
