import sys
from pathlib import Path

sys.path.append("/Users/hgz/Documents/researchHBLab/Projects/Selex-X-Genome/source")

from base import DiskFile, ENCODE_Object, Experiment
from disk_files import PE_Fastq, SE_Fastq
from experiment import Control, TFChipSeq
from files import PE_File, SE_File
from library import Library
from search import EncodeSearch

Human_MYC = EncodeSearch("MYC", "Homo+sapiens", limit="1")
Human_MYC.search()

# print(Human_MYC.search_result)
Experiments = Human_MYC.get_experiments()

Experiments[0].fetchData()
libraries = Experiments[0].get_libraries()

files = libraries[0].get_Files()

print(files[0].href)
files[0].download(Path("/Users/hgz/Documents/researchHBLab/Projects/Selex-X-Genome"))
