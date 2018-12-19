from math import sqrt
from progressbar import ProgressBar

def dist(spec1, spec2):
    """
    Compute distance between two spectra.

    Parameters
    ----------
    spec1, spec2: dict
        Spectra.

    Output
    ------
    Distance
    """
    k1 = set(spec1.keys())
    k2 = set(spec2.keys())
    d = 0
    for k in k1:
        if k in k2:
            d += sqrt((spec1[k] - spec2[k])**2)
        else:
            d += spec1[k]
    for k in (k2 - k1):
        d += spec2[k]
    return d

def running_dist(k, seq, full_spec, win):
    L = len(seq)
    win_specs = [dict() for i in range(win)]
    bar = ProgressBar()
    distances = list()
    for i in bar(range(L-k)):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            if i - win < 0:
                for spec in win_specs[:i+1]:
                    if kmer not in spec:
                        spec[kmer] = 1/win
                    else:
                        spec[kmer] += 1/win
            else:
                distances.append(dist(full_spec, win_specs.pop(0)))
                win_specs.append(dict())
                for spec in win_specs[i-win+1:i+1]:
                    if kmer not in spec:
                        spec[kmer] = 1/win
                    else:
                        spec[kmer] += 1/win

    return distances
