import copy
from multiprocessing import Process, Manager

  #################
## SPECTRUM OBJECT ############################################################
  #################

class Spectrum:
    """
    K-mer spectrum container.
    """

    # INITIALIZATION ----------------------------------------------------------
    def __init__(self, kmer_freqs=None, normalized=False):
        """
        Constructor

        Parameters
        ----------
        kmer_freqs: dict
            Dictionary of K-mer frequencies built as:
                {<kmer_str>: <frequency>}

        normalized: boolean
            Set it to True if kmer_freqs is normalized.
        """
        if kmer_freqs is None:
            self.kmer_freqs = dict()
        else:
            self.kmer_freqs = kmer_freqs
        self.normalized = normalized


    # TOOLS -------------------------------------------------------------------
    def normalize(self, coef):
        for k in self.kmer_freqs:
            self.kmer_freqs[k] = self.kmer_freqs[k] / coef
        self.normalized = True

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

    def __add__(self, other):
        if type(other) != Spectrum:
            err_str = "Addition not implemented yet for %s" % type(other)
            raise Exception(err_str)
        if self.normalized or other.normalized:
            print("WARNING: you are adding normalized spectra !")
        new_spectrum = copy.deepcopy(self)
        for k in other.kmers():
            if k in self:
                new_spectrum[k] = self[k] + other[k]
            else:
                new_spectrum[k] = other[k]
        return new_spectrum


  #########################
## SPECTRUM CREATION TOOLS ####################################################
  #########################

def kmer_spectrum(k, seq, normalize=True):
    """
    Get the k-mer frequencies.

    Parameters
    ----------
    k: int
        Number of amino-acids in K-mers.

    full_seq: string
        DNA sequence.
    
    normalize: boolean
        If True normalize by the length of the sequence.

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
    if normalize:
        spec.normalize(L-k+1)

    return spec

def mproc_kmer_spectrum(k, seq, normalize=True, n_process=4):
    """
    Get the k-mer frequencies.

    Parameters
    ----------
    k: int
        Number of amino-acids in K-mers.

    full_seq: string
        DNA sequence.

    normalize: boolean
        If True normalize by the length of the sequence.

    Output
    ------
    Normalized spectrum object
    """
    # Initialization ----------------------------------------------------------
    spec = Spectrum()
    L = len(seq)
    L_process = (L-k) // n_process
    jobs = list()
    manager = Manager()
    results = manager.list([None]*n_process)

    # Iterate on tasks --------------------------------------------------------
    for i in range(n_process):

        # > If last process
        if i == n_process - 1:
            stop = L
        # > For the other processes
        else:
            stop = (i + 1)*L_process + k - 1

        start = i*L_process

        # > Get sequences for processes
        trunc_seq = seq[start:stop]

        # > Launch task
        job = Process(target=worker,
                      args=(results, i, k, trunc_seq, False))
        job.start()
        jobs.append(job)

    # Iterate on tasks --------------------------------------------------------
    for job in jobs:
        # > Wait the end of the execution
        job.join()

    # > Retrieve results
    for r in results:
        spec += r

    if normalize:
        spec.normalize(L-k+1)

    return spec

def worker(results, i, *args, **kwargs):
    results[i] = kmer_spectrum(*args, **kwargs)

