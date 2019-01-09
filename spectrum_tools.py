from math import sqrt
from spectrum import Spectrum
from progressbar import ProgressBar
from threading import Thread

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

    # Get K-mers --------------------------------------------------------------
    k1 = set(spec1.keys())
    k2 = set(spec2.keys())
    
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

def running_dist(k, seq, full_spec, win, step=1):
    """
    Compute distance on a running window.

    Parameters
    ----------
    k: int
        Number of amino acid in the K-mers.

    seq: str
        Sequence to use.

    full_spec: spectra
        Spectra of the full sequence.

    win: int
        Window length.

    step: int
        Step between windows.

    Output
    ------
    Distances as a list of floats.
    """
    # Initialization ----------------------------------------------------------
    L = len(seq)
    win_specs = [Spectrum() for i in range(win//step)]
    bar = ProgressBar()
    distances = list()

    # Iterate on amino K-mers of seq ------------------------------------------
    for i in bar(range(L-k)):
        kmer = seq[i:i+k]
        # > Discard K-mer containing not defined amino acids
        if "N" not in kmer:
            # > Left border
            if i - win < 0:
                # Iterate on windows containing the current K-mer -------------
                win_i = i // step
                for idx, spec in enumerate(win_specs[:win_i+1]):
                    if kmer not in spec:
                        spec[kmer] = 1/win
                    else:
                        spec[kmer] += 1/win
            else:
                # > If a window is full, drop it and compute distance
                if (i - win) % step == 0:
                    win_spec = win_specs.pop(0)
                    distances.append(dist(full_spec, win_spec))
                    win_specs.append(Spectrum())

                # Iterate on windows containing the current K-mer -------------
                for spec in win_specs:
                    if kmer not in spec:
                        spec[kmer] = 1/win
                    else:
                        spec[kmer] += 1/win
    return distances

def multithread_running_dist(k, seq, full_spec, win, step=1, n_thread=4):
    L = len(seq)
    N = (L - k) // n_thread
    thds = list()
    distances = list()
    for i in range(n_thread):
        if i == n_thread - 1:
            stop = L - k
        else:
            stop = (i + 1)*N
        thd = RunningDistThread(i*N, stop, k, seq, full_spec, win, step)
        thd.start()
        thds.append(thd)
    for thd in thds:
        thd.join()
        if thd.distances == None:
            raise Exception("FUUUUCK OF")
        distances += thd.distances

    return distances

class RunningDistThread(Thread):
    def __init__(self, start, stop, k, seq, full_spec, win, step):
        Thread.__init__(self)
        self.start_ = start
        self.stop = stop
        self.k = k
        self.seq = seq
        self.full_spec = full_spec
        self.win = win
        self.step = step
        self.distances = None

    def run(self):
        L = len(self.seq)
        win_specs = [dict() for i in range(self.win//self.step)]
        bar = ProgressBar()
        distances = list()
        for i in range(self.start_, self.stop + self.win):
            kmer = self.seq[i:i+self.k]
            if "N" not in kmer:
                if i - self.win < 0:
                    for spec in win_specs[:i+1]:
                        if kmer not in spec:
                            spec[kmer] = 1/self.win
                        else:
                            spec[kmer] += 1/self.win
                else:
                    if (i - self.win) % self.step == 0:
                        distances.append(dist(self.full_spec, win_specs.pop(0)))
                        win_specs.append(dict())
                    for spec in win_specs:
                        if kmer not in spec:
                            spec[kmer] = 1/self.win
                        else:
                            spec[kmer] += 1/self.win
        self.distances = distances

def running_average(X, win):
    cumsum = list()
    rslt = list()
    for i, x in enumerate(X):
        if i == 0:
            cumsum.append(x)
        else:
            cumsum.append(cumsum[-1] + x)
        if i + 1 > win:
            rslt.append((cumsum[-1] - cumsum[-win])/win)
    return rslt
