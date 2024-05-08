class Motif:
    def __init__(self, tf: str, organism: str):
        self.tf = tf
        self.organism = organism

    def score_seq(self):
        """Score a sequence"""
        pass

    def score_random_kmers(self):
        """Get the average score and score distribution of the motif scored against data/random_sequences.csv"""
        pass
