from ENCODE.search import EncodeSearch
from motifs.motif import Mononucleotide
from motifs.parse_motifcentral_json import MOTIFCENTRAL

motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[867])

search_result = EncodeSearch.search_using_MOTIFCENTRAL_json(MOTIFCENTRAL[867])

motif.motif_X_TF(search_result)


# ! I need a way to ensure that each job ran successfully.
