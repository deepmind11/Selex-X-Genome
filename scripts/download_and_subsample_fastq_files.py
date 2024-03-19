#!/usr/bin/env python3

import json
import random
import subprocess
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests


def process_fastqgz_file(fastq_file_path, size=2 * (10**6)):
    """
    0.1 Checks if the seed_log_file exists, if it does reads seed from seed_log_file
    0. else , sets random seed for subsampling and saves it to a seed_log_file

    1. Subsamples a fastq.gz file and deletes the original
    2. Converts the subsample file to fasta format and deletes the original

    Args:
        fastq_file_path (Path): Path to the .fastq.gz file
        size (int, optional): Size of the sample. Defaults to 2 * (10**6).
    """

    # Defining the seed path
    seed_Path = fastq_file_path.parent / Path("seed.txt")

    # Read/Set the seed
    if seed_Path.exists():
        with seed_Path.open(mode="r") as seed_file:
            seed = int(seed_file.read_text())
    else:
        with seed_Path.open(mode="w") as seed_file:
            seed = random.randint(1, 100)
            seed_file.write(str(seed))

    # Defining Path to subsampled fastq and fasta files
    subsampled_fastq_file = fastq_file_path.parent / Path(
        fastq_file_path.stem + "_subsampled.fq"
    )
    subsampled_fasta_file = fastq_file_path.parent / Path(
        fastq_file_path.stem + "_subsampled.fa"
    )

    # Calling seqtk to subsample and convert to fasta
    subprocess.run(
        [
            "seqtk",
            "sample",
            f"-s{seed}",
            str(fastq_file_path),
            str(size),
            ">",
            str(subsampled_fastq_file),
        ]
    )
    subprocess.run(
        [
            "seqtk",
            "seq",
            "-a",
            str(subsampled_fastq_file),
            ">",
            str(subsampled_fasta_file),
        ]
    )

    # Deleting the original and subsampled fastq files
    fastq_file_path.unlink()
    subsampled_fastq_file.unlink()

    return


#
def download_pe_files(library, files, experiment_json_DR):
    """Downloads fastq files (one pair for PE) for a given library (belonging to a particular experiment).

    Args:
        library (string): ENCODE accession of the library ENCLBXXXXXX
        files (list):     A list of fastq file objects(jsons) belonging to the library.
        experiment_json_DR (Path): Path to the directory containing the experiment json file.
    """

    # Choosing the files
    file_to_download = random.choice(files)
    ## Getting the corresponding paired end file
    for file in files:
        if file["accession"] == file_to_download["paired_with"][7:][:-1]:
            pe_file = file
            break

    # Getting the download url
    download_url = "https://www.encodeproject.org" + file_to_download["href"]
    download_url_pe = "https://www.encodeproject.org" + pe_file["href"]

    # Constructing the filnames and filepaths
    filename = (
        f'{file_to_download["accession"]}_{file_to_download["paired_end"]}.fastq.gz'
    )
    filepath = experiment_json_DR / Path(library) / Path(filename)
    ## For PE
    filename_pe = f'{pe_file["accession"]}_{pe_file["paired_end"]}.fastq.gz'
    filepath_pe = experiment_json_DR / Path(library) / Path(filename_pe)

    filepath.parent.mkdir(parents=True, exist_ok=True)
    filepath_pe.parent.mkdir(parents=True, exist_ok=True)

    # Downloading the PE files in chunks
    response = requests.get(download_url, stream=True)
    with filepath.open(mode="wb") as file:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            file.write(chunk)
    response_pe = requests.get(download_url_pe, stream=True)
    with filepath_pe.open(mode="wb") as file_pe:
        for chunk_pe in response_pe.iter_content(chunk_size=1024 * 1024):
            file.write(chunk_pe)

    return


def download_se_files(library, files, experiment_json_DR):
    """Downloads fastq files (one for SE) for a given library (belonging to a particular experiment).

    Args:
        library (string): ENCODE accession of the library ENCLBXXXXXX
        files (list):     A list of fastq file objects(jsons) belonging to the library.
        experiment_json_DR (Path): Path to the directory containing the experiment json file.
    """

    # Choosing the file
    file_to_download = random.choice(files)

    # Getting the download url
    download_url = "https://www.encodeproject.org" + file_to_download["href"]

    # Constructing the filnames and filepaths
    filename = f'{file_to_download["accession"]}.fastq.gz'
    filepath = experiment_json_DR / Path(library) / Path(filename)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # Downloading the file in chunks
    response = requests.get(download_url, stream=True)
    with filepath.open(mode="wb") as file:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            file.write(chunk)

    return


def download_and_subsample_fastq_files(experiment_json):
    """
    1. Download fastq files for all libraries belong to an experiment.
       Downlaod location: ENCSRXXXXXX/ENCLBXXXXXX/file.fastq.gz
    2. Subsample the fastq.gz files.
    3. Compress the subsampled file.
    4. Delete the original downloaded file.

    Args:
        experiment_json (Path): Path to the experiment_data json file.
    """
    # Loading the json data in python dict.
    with open(experiment_json, "r") as f:
        expr_data = json.load(f)

    libraries = list(
        [
            expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(expr_data["replicates"]))
        ]
    )

    # Gather fastq files for each library
    library_fastq_files = dict(zip(libraries, [[] for i in range(len(libraries))]))

    for file in expr_data["files"]:
        if "replicate" in file.keys() and file["file_format"] == "fastq.gz":
            library_fastq_files[file["replicate"]["library"][11:-1]].append(file)

    # Download the files for each library concurrently

    # Argument list to passed to executor.map
    arguments = list(
        [
            (library, files, Path(experiment_json).parent)
            for library, files in library_fastq_files.items()
        ]
    )

    if (
        library_fastq_files[list(library_fastq_files.keys())[0]][0]["run_type"]
        == "paired-ended"
    ):
        with ThreadPoolExecutor() as executor:
            executor.map(lambda args: download_pe_files(*args), arguments)
    else:
        with ThreadPoolExecutor() as executor:
            executor.map(lambda args: download_se_files(*args), arguments)

    # Process the fastq files (Using seqtk)
    for library in library_fastq_files:

        # Directory of the fastq files
        fastq_files_dir = Path(experiment_json).parent / Path(str(library))

        # Process if fastq file
        for file in fastq_files_dir.iterdir():
            if file.suffix == ".gz":
                process_fastqgz_file(file)

    return


if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("experiment_json", type=str)
    args = parser.parse_args()

    download_and_subsample_fastq_files(Path(args.experiment_json))
