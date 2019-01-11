from threading import Thread

from progressbar import ProgressBar

from .tools import dist
from .spectrum import Spectrum

  ################
## DISTANCE TOOLS #############################################################
  ################

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
    win_specs = []
    distances = list()

    # Iterate on K-mers of seq ------------------------------------------------
    for i in range(L-k+1):
        kmer = seq[i:i+k]

        # > Discard K-mers containing not defined amino acids
        if "N" not in kmer:

            # > Add window each steps
            if i % step == 0:
                win_specs.append(Spectrum())

            # Iterate on windows containing the current K-mer -----------------
            for spec in win_specs:
                if kmer not in spec:
                    spec[kmer] = 1.
                else:
                    spec[kmer] += 1.

            # > If a window is full, drop it and compute distance
            if (i - win + 1) % step == 0 and i >= win-1:
                win_spec = win_specs.pop(0)
                win_spec.normalize(win)
                distances.append(dist(full_spec, win_spec))

    return distances

  ############################
## MULTITHREAD DISTANCE TOOLS #################################################
  ############################

def multithread_running_dist(k, seq, full_spec, win, step=1, n_thread=4):
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

    n_thread: int
        Number of thread to use

    Output
    ------
    Distances as a list of floats.
    """
    # Initialization ----------------------------------------------------------
    L = len(seq)
    L_thread = (L-k-win) // n_thread
    thds = list()
    distances = list()

    # Iterate on threads ------------------------------------------------------
    for i in range(n_thread):

        # > If last thread
        if i == n_thread - 1:
            stop = L
        # > For the other threads
        else:
            stop = (i + 1)*L_thread + k + win - 2

        if i == 0:
            start = 0
        else:
            if i*L_thread % step != 0:
                start = i*L_thread + step - ((i*L_thread) % step)
            else:
                start = i*L_thread

        # > Get sequences for threads
        trunc_seq = seq[start:stop]

        # > Launch thread
        thd = RunningDistThread(k, trunc_seq, full_spec, win, step)
        thd.start()
        thds.append(thd)

    # Iterate on treads -------------------------------------------------------
    for thd in thds:
        thd.join()
        if thd.distances == None:
            raise Exception("Error in one thread")

        # > Retrieve distances
        distances += thd.distances

    return distances

class RunningDistThread(Thread):
    """
    Running distance thread.
    """
    def __init__(self, *args, **kwargs):
        Thread.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.distances = None

    def run(self):
        self.distances = running_dist(*self.args, **self.kwargs)

  #########
## FILTERS ####################################################################
  #########

def running_average(X, win):
    """
    Low pass filter.

    Parameters
    ----------
    X: Array-like
        Signal to filter

    win: int
        Window length. Bigger is the window, lower is the cutoff frequency.

    Output
    ------
    Filtered signal.
    """
    cumsum = list()
    rslt = list()
    for i, x in enumerate(X):
        if i == 0:
            cumsum.append(x)
        else:
            cumsum.append(cumsum[-1] + x)
        if i + 1 > win:
            rslt.append((cumsum[-1] - cumsum[-win -1]) / win)
    return rslt
