#!/usr/bin/env python3

"""
Given a TF and organism. Downloads and transforms all the fastq files from ENCODE. 
Transfrom it to fasta based on source/se_processing_template.sh
The script is step 1 in the pipeline it is followed by the creation of the count tables.
"""


from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from ENCODE.experiment import Control, TFChipSeq
from ENCODE.library import Library
from ENCODE.search import EncodeSearch


def process_library(library: Library, libpath: Path):
        """
        Download all the fastq files for a library. Process it and delete the originals.
        
        library: ENCODE library object.
        libpath: Path to the basedr for given library.

        """
        
        # All fastq files associated with the library
        files = library.get_Files()

        for file in files:
            # if fasta file exists then skip loop (checks if gz fasta exists)
            if (libpath / Path (f'{file.accession}.fasta')).exists():
                continue
            else:
                # dowload the fastq files
                fq = file.download(libpath)
                # Subsample the fastq file. Default is half a mil
                sub_fq = fq.subsample()
                # Delete the og file
                fq.delete()
                # Submit job on slurm to transfrom file
                sub_fq.slurm_transform()


def process_controls(controls: list[Control]):
    """Downloads all the fastq files for the given set of control experiments. Transform the fastq files to fasta.
        Delete all fastq files.

    Args:
        controls (list): List of encode Control Objects.
    """

    for control in controls:
        try:
            # fetch data from ENCODE
            control.fetchData()
            # Get the libraries
            libraries = control.get_libraries()

            # If control already exists, skip iteration
            if Path(f"/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/Control/{control.accession}").exists():
                continue
            else:
                # Library is uniquely identified by a technical and biological replicate
                for library in libraries:
                    process_library(library, Path(f'/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/Control/{control.accession}/{library.accession}'))
                    
        except Exception as e:
            print(f"Error: {control.accession}")
            print(e)


def process_exp(experiment: TFChipSeq):
    """For a given TFChipSeq object. It download all the important files along with the controls.
        It then process the fastq files. And calls a slurm job to transform it to fasta.
        Delete all fastq files.

    Args:
        experiment (TFChipSeq): A Chip seq ENCODE experiment.
    """
    try:
        
        # Fetch data from ENCODE
        experiment.fetchData()
        
        # Get libraries and controls
        libraries = experiment.get_libraries()
        controls = experiment.get_controls()
        
        # Process the controls
        process_controls(controls)

        # Library is uniquely identified by a technical and biological replicate #
        for library in libraries:
            
            libpath = Path(f"/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/{experiment.accession}/{library.accession}")

            # If library exists skip iteration
            if libpath.exists():
                continue
            else:
                process_library(library, libpath)
                    
          

    except Exception as e:
        print(f"Error: {experiment.accession}")
        print(e)




if __name__ == "__main__":

    # Parsing the arguments
    parser = ArgumentParser()
    parser.add_argument("TF", type=str)
    parser.add_argument("organism", type=str)  # Homo+sapiens
    args = parser.parse_args()

    # Check that the organism has been specified in the right format.
    if "+" not in args.organism:
        raise ValueError("Organism must be in the format of Homo+sapiens")

    # Querying the encode DB
    query = EncodeSearch(args.TF, args.organism, limit="all")
    query.search()

    # All the experiments for the query
    Experiments = query.get_experiments()

    # Processing the experiment. All data is downloaded to data/TF_Homo_sapiens
    with ThreadPoolExecutor() as executor:
        executor.map(process_exp, Experiments)
