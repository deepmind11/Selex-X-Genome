from argparse import ArgumentParser
from pathlib import Path

from utils.slurmjob import Slurmjob

if __name__ == "__main__":

    # Parsing the arguments
    parser = ArgumentParser()
    parser.add_argument("TF", type=str)
    parser.add_argument("organism", type=str)  # Homo+sapiens
    parser.add_argument("motifcentral_index", type=int)
    args = parser.parse_args()

    # Downloading the files from ENCODE
    slurmjob_path_1 = Path(
        f'/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/download.job'
    )
    job_name_1 = "Download_files_from_ENCODE"
    job_script_1 = Path(
        "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/source/analysis/tf_organism_1.py"
    )
    job_params_1 = (args.TF, args.organism)
    job_cores_1 = 6
    job_time_1 = 6

    # Create slurm job file
    job_file_1 = Slurmjob(
        file_path=slurmjob_path_1,
        job_name=job_name_1,
        job_script=job_script_1,
        job_params=job_params_1,
        cores=job_cores_1,
        time=job_time_1,
    )
    job_file_1.create_file()
    job_file_1.submitJob()

    # Scoring the count tables
    slurmjob_path_2 = Path(
        f'/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/{args.TF}_{args.organism.replace("+","_")}/{args.motifcentral_index}_score.job'
    )
    job_name_2 = "Score_count_tables"
    job_script_2 = Path(
        "/burg/hblab/users/hg2604/Projects/Selex-X-Genome/source/analysis/tf_organism_2.py"
    )
    job_params_2 = (args.TF, args.organism, args.motifcentral_index)
    job_cores_2 = 32
    job_time_2 = 11

    # Create slurm job file
    job_file_2 = Slurmjob(
        file_path=slurmjob_path_2,
        job_name=job_name_2,
        job_script=job_script_2,
        job_params=job_params_2,
        cores=job_cores_2,
        time=job_time_2,
    )
    job_file_2.create_file()
    # job_file_2.submitJob()
