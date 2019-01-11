from math import sqrt
from .spectrum import Spectrum
from progressbar import ProgressBar

def load_file(path, verbose=False):
    """
    Load fna file.

    Parameter
    ---------
    path: string
        Path to the file
    
    verbose: boolean
        If True print some informations on the sequence.

    Output
    ------
    Sequence as a string.
    """
    seq = str()
    with open(path, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i == 0:
                continue
            seq += line
    seq = seq.replace('\n', '')
    if verbose:
        print('Length of the genome:', len(seq))
    return seq

def dist(spec1, spec2):
    """
    Compute distance between two spectra.

    Parameters
    ----------
    spec1, spec2: Spectrum objects

    Output
    ------
    Distance
    """

    # Get K-mers --------------------------------------------------------------
    k1 = set(spec1.kmers())
    k2 = set(spec2.kmers())
    
    # > d = distance
    d = 0

    # Iterate on K-mers of spec1 ----------------------------------------------
    for k in k1:
        if k in k2:
            d += sqrt((spec1[k] - spec2[k])**2)
        else:
            d += spec1[k]

    # Iterate on K-mers of spec2 that are not in spec1 ------------------------
    for k in (k2 - k1):
        d += spec2[k]

    return d

