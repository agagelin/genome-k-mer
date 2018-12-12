from math import sqrt

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
