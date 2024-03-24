from pathlib import Path
import subprocess


class Fastq:
    def __init__(self, file_path: Path = None):
        self.file_path = file_path

    def delete(self):
        """Deletes the fastq file from disk."""
        self.file_path.unlink()
        self.file_path = None
    
    def subsample(self, seed, size):
        """Returns subsampled fastq object"""
        
        subsampled_fastq_file = self.file_path.parent / Path(self.file_path.name[:11] + "_subsampled.fq")
        
        with open(subsampled_fastq_file, "w") as output_file:
            subprocess.run(
                ["seqtk", "sample", f"-s{seed}", str(self.file_path), str(size)],
                stdout=output_file,
            )
        
        return Fastq(subsampled_fastq_file)
    
    def transform(self, transform_script):

        


