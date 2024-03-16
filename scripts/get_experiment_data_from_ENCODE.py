"""
This program queries the ENCODE DB using its REST API. Takes as input the search results.

1. Saves all the data for an experiment in a json. data/tf_organism/experiments/ENCSRXXXXXX/experiment_data.json
2. Saves data for the control experiment as a json. data/control/ENCSRXXXXXX/exepriment_data.json

3. Saves only one (picks one randomly, if multiple available. Downloads 1 pair for PE) fastq files belong to a particular 
   library (defined by tech replicate # & bio replicate #) at
   data/tf_organism/experiments/ENCSRXXXXXX/ENCLBXXXXXX/file.fastq

   3.1) Subsamples the fastq file 2,000,000. gzips it
   3.2) deletes the downloaded file.
"""

import json
import random
from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path
from urllib.request import urlretrieve

import requests


# helper functions
def download_files(library_fastq_files,experiment_json_DR):
    """Downloads fastq files (one for SE, one pair for PE) for a given experiment.

    Args:
        library_fastq_files (dict): A dictionary of all the libraries in an experiment. Each library has a list of fastq file objects.
        experiment_json_DR (Path): Path to the directory containing the experiment json file.
    """
    if library_fastq_files[list(library_fastq_files.keys())[0]]['run_type'] == 'paired-ended':

        # Download the PE files
        for library in library_fastq_files.keys():
            
            # Choosing the files
            file_to_download = random.choice(library_fastq_files[library])
            ## Getting the corresponding paired end file
            for file in library_fastq_files[library]:
                if file['accession'] == file_to_download['paired_with'][7:][:-1]:
                    pe_file = file
                    break
            
            # Getting the download url
            download_url = "https://www.encodeproject.org" + file_to_download['href']
            download_url_pe = "https://www.encodeproject.org" + pe_file['href']
            
            # Constructing the filnames and filepaths
            filename = f'{file_to_download['accession']}_{file_to_download['paired_end']}.fastq'
            filepath = experiment_json_DR / Path(library) / Path(filename)
            ## For PE
            filename_pe = f'{pe_file['accession']}_{pe_file['paired_end']}.fastq'
            filepath_pe = experiment_json_DR / Path(library) / Path(filename_pe)

            filepath.parent.mkdir(parents=True, exist_ok=True)
            filepath_pe.parent.mkdir(parents=True, exist_ok=True)
          
            # Downloading the PE files in chunks
            response = requests.get(download_url, stream=True)
            with filepath.open(mode = 'wb') as file:
                for chunk in response.iter_content(chunk_size = 1024*1024):
                    file.write(chunk)
            response_pe = requests.get(download_url_pe, stream=True)
            with filepath_pe.open(mode = 'wb') as file_pe:
                for chunk_pe in response_pe.iter_content(chunk_size = 1024*1024):
                    file.write(chunk_pe)
      
    
    else:
        # Downloading the SE files
        for library in library_fastq_files.keys():
            file_to_download = random.choice(library_fastq_files[library])
            download_url = "https://www.encodeproject.org" + file_to_download['href']
            filename = f'{file_to_download['accession']}.fastq'
            filepath = experiment_json_DR / Path(library) /Path(filename)
            filepath.parent.mkdir(parents=True, exist_ok=True)
            
            #Downloading the file in chunks
            response = requests.get(download_url, stream=True)
            with filepath.open(mode = 'wb') as file:
                for chunk in response.iter_content(chunk_size = 1024*1024):
                    file.write(chunk)
            
        
    return


def save_experiment_data_as_json(search_result, base_dir):
    """
    Queries the ENCODE DB using the experiment accession and saves the json data at {target_dir}.

    Args:
        search_result (dict: @graph): An element of the @graph list returned by ENCODE search
        base_dir (Path):       The base directory where all the experiments will be saved.
    """

    expr_accession = search_result["accession"]

    # Force return from the server in JSON format
    headers = {"accept": "application/json"}

    # This searches the ENCODE database for the phrase "bone chip"
    url = f"https://www.encodeproject.org/experiments/{expr_accession}"

    # GET the search result
    response = requests.get(url, headers=headers)

    # Extract the JSON response as a python dictionary
    expr_data = response.json()

    # Constructing the path to save the json
    target_dir = base_dir / Path(expr_accession)
    target_dir.mkdir(parents=True, exist_ok=True)

    # Saving search results as json
    json_file = target_dir / Path("experiment_data.json")
    json_file.touch()

    with json_file.open(mode="w") as js:
        json.dump(expr_data, js, indent=4)

    # Logging Info. Storing all the libraries and the control belonging to an experiment
    log_file = target_dir / Path("experiment.log")
    log_file.touch()

    libraries = list(
        [
            expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(expr_data["replicates"]))
        ]
    )
    controls = list(
        [expr_data["possible_controls"][i]["accession"]]
        for i in range(len(expr_data["possible_controls"]))
    )

    with log_file.open(mode="r") as log:
        
        log.write(f"The query paramerters are tf={args.tf} organism={args.organism} expr_accession={expr_accession}{chr(10)}")
        log.write(f"Date: {datetime.now():%c}{chr(10)}")
        log.write(f"List of Libraries: {*libraries,}{chr(10)}")
        log.write(f"Possible Controls: {*controls,}{chr(10)}")

    return 



# ! Need to write code to subsample the fastq files.
# ! Also, need a way to parallelize downloads.
def download_and_subsample_fastq_files(experiment_json):
    """
    1. Download fastq files for all libraries belong to an experiment.
       Downlaod location: ENCSRXXXXXX/ENCLBXXXXXX/file.fastq.gz
    2. Subsample the fastq.gz files. Decompress and subsample. 2 million reads.
    3. Compress the subsampled file.
    4. Delete the original downloaded file.

    Args:
        experiment_json (file(.json)): Path to the experiment_data json file.
    """
    #Loading the json data in python dict.
    with open(experiment_json, "r") as f:
        expr_data = json.load(f)


    libraries = list(
        [
            expr_data["replicates"][i]["library"]["accession"]
            for i in range(len(expr_data["replicates"]))
        ]
    )

    #Gather fastq files for each library
    library_fastq_files = dict(zip(libraries,[[] for i in range(len(libraries))]))

    for file in  expr_data['files']:
        if 'replicate' in file.keys() and file['file_format']=='fastq':
            library_fastq_files[file['replicate']['library'][11:-1]].append(file)

    # Download the files
    download_files(library_fastq_files,Path(experiment_json).parent)

    # Subsample the fastq files (Using seqtk)

    # Compress the subsampled file
    # Delete the original downloaded file

    return

        
def save_control_data_for_experiment(experiment_json):
    """
    Takes as input an exepriment_data.json file. Check if the control experiment data exists, if not 
    download it. And then proceed to download all the relevant files for that control experiment.

    Args:
        experiment_json (file(.json)): data/tf_organism/experiments/ENCSRXXXXXX/experiment_data.json
    """



if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("tf", type=str)
    parser.add_argument("organism", type=str)

    args = parser.parse_args()

    # loading the search results
    search_result_json = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}/search_results.json'
    )

    with search_result_json.open(mode="r") as srj:
        search_results = json.load(srj)

    search_results = search_results["@graph"]

    # Querying the ENCODE DB for all the hits(expriments) and saving as json.
    base_dir = Path(__file__).parent.parent / Path(
        f'data/{args.tf}_{args.organism.replace("+", "_")}/experiments'
    )

    # ! I should implement concurrency here. 
    for search_result in search_results:
        save_experiment_data_as_json(search_result, base_dir)

        # !Downlaod the files for that experiment

        # !Downlaod the control files for that experiment


