  #################
## SPECTRUM OBJECT ############################################################
  #################

class Spectrum:
    """
    K-mer spectrum container.
    """

    # INITIALIZATION ----------------------------------------------------------
    def __init__(self, kmer_freqs=None):
        """
        Constructor

        Parameters
        ----------
        kmer_freqs: dict
            Dictionary of K-mer frequencies built as:
                {<kmer_str>: <frequency>}
        """
        if kmer_freqs is None:
            self.kmer_freqs = dict()
        else:
            self.kmer_freqs = kmer_freqs

    # TOOLS -------------------------------------------------------------------
    def normalize(self, coef):
        for k in self.kmer_freqs:
            self.kmer_freqs[k] = self.kmer_freqs[k] / coef

    # ACCESS DATA -------------------------------------------------------------
    def kmers(self):
        return list(self.kmer_freqs.keys())

    def kmer_number(self):
        return len(self.kmers())

    # PRINT OVERLOAD ----------------------------------------------------------
    def __str__(self):
        kmers = self.kmers()
        L = len(kmers)
        K = len(kmers[0])
        str_ = "<Spectrum: {} {}-mers>".format(L, K)
        return str_

    # OPERATOR OVERLOAD ------------------------------------------------------- 
    def __getitem__(self, key):
        return self.kmer_freqs[key]

    def __setitem__(self, key, items):
        self.kmer_freqs[key] = items

    def __contains__(self, key):
        return key in self.kmer_freqs


  #########################
## SPECTRUM CREATION TOOLS ####################################################
  #########################

def get_kmer_spectrum(k, seq):
    """
    Get the k-mer frequencies.

    Parameters
    ----------
    k: int
        Number of amino-acids in K-mers.

    full_seq: string
        DNA sequence.

    Output
    ------
    Spectrum object
    """
    spec = Spectrum()
    L = len(seq)
    for i in range(L-k+1):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            if kmer not in spec:
                spec[kmer] = 1.
            else:
                spec[kmer] += 1.
    spec.normalize(L-k+1)

    return spec

