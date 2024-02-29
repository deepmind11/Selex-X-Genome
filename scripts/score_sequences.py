import numpy as np

def score_window(numpy_motif, seq):
    """Scores a window the size of motif 

    Args:
        numpy_motif (numpy array of shape len(seq) x 4): The motif
        seq (String): The DNA Sequence "ACGT" szie of the motif
    Returns:
        score (float): The score of the window
    """
    score = 0
    for i, base in enumerate(seq):
        match base:
            case 'A':
                score += numpy_motif[0, i]
            case 'C':
                score += numpy_motif[1, i]
            case 'G':
                score += numpy_motif[2, i]
            case 'T':
                score += numpy_motif[3, i]
    
    return score

def score_seq_against_motif(seq, motif):
    """
    Given a sequence and a motif, scores the sequence against the motif.

    arguments:
        1. seq (str): sequence to score (ACGT)
        2. motif (list): A list of deldelGs. Lenght = Size of motif * 4. Order = ACGT
    returns:
        1. score (float): the score of the sequence against the motif in log space.
    """
    

    motif_size = int(len(motif)/4)
    numpy_motif = np.array(motif).reshape(4, motif_size, order='F')

    score = 0
    for i in range(0, len(seq)-2):
        score += score_window(numpy_motif, seq[i:i+motif_size])
    
    return score
