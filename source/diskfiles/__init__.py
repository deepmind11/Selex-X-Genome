__all__ = [
    "DiskFile",
    "CountTable",
    "PE_Fastq",
    "SE_Fastq",
    "Slurmjob",
]


from .base import DiskFile
from .countTables import CountTable
from .fastq import PE_Fastq, SE_Fastq
from .slurmjob import Slurmjob

# del base, countTables, fastq, slurmjob
