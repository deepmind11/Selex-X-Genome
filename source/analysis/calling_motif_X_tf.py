from ENCODE.search import EncodeSearch
from motifs.motif import Mononucleotide
from motifs.parse_motifcentral_json import MOTIFCENTRAL

motif = Mononucleotide.create_from_motif_central(MOTIFCENTRAL[465])

search_result = EncodeSearch.search_using_MOTIFCENTRAL_json(MOTIFCENTRAL[628])
# search_result = EncodeSearch('MAX', 'Homo+sapiens')


motif.motif_X_TF(
    search_result, data_path="/burg/hblab/users/wl2924/Selex-X-Genome/data/"
)

# ! I need a way to ensure that each job ran successfully.
