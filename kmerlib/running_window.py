from threading import Thread
from multiprocessing import Process, Manager
from time import sleep

from .tools import *
from .spectrum import Spectrum
from utils.term_colors import *

  ################
## DISTANCE TOOLS #############################################################
  ################

def running_dist(k, seq, full_spec, win, step=1, fast=True):
    """
    Compute distance on a running window.

    Parameters
    ----------
    k: int
        Number of nucleic acid in the K-mers.

    seq: str
        Sequence to use.

    full_spec: spectra
        Spectra of the full sequence.

    win: int
        Window length.

    step: int
        Step between windows.

    fast: boolean
        Use fast mode.

    Output ------
    Distances as a list of floats.
    """
    # Initialization ----------------------------------------------------------
    L = len(seq)
    win_specs = []
    distances = list()

    # Iterate on K-mers of seq ------------------------------------------------
    for i in range(L-k+1):
        kmer = seq[i:i+k]

        # > Discard K-mers containing not defined nucleic acids
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
                if fast:
                    distances.append(faster_dist(full_spec, win_spec))
                else:
                    distances.append(dist(full_spec, win_spec))

    return distances


  #############################
## MULTIPROCESS DISTANCE TOOLS ################################################
  #############################

def mproc_running_dist(k, seq, full_spec, win, step=1, fast=True, n_process=4):
    """
    Compute distance on a running window.

    Parameters
    ----------
    k: int
        Number of nucleic acid in the K-mers.

    seq: str
        Sequence to use.

    full_spec: spectra
        Spectra of the full sequence.

    win: int
        Window length.

    step: int
        Step between windows.

    fast: boolean
        Use fast mode.

    n_process: int
        Number of process to use

    Output
    ------
    Distances as a list of floats.
    """
    # Initialization ----------------------------------------------------------
    L = len(seq)
    L_process = (L-k-win) // n_process
    jobs = list()
    distances = list()
    manager = Manager()
    results = manager.list([None]*n_process)

    # Iterate on tasks --------------------------------------------------------
    for i in range(n_process):

        # > If last process
        if i == n_process - 1:
            stop = L
        # > For the other processes
        else:
            stop = (i + 1)*L_process + k + win - 2

        # > If first process
        if i == 0:
            start = 0
        # > For the others
        else:
            if i*L_process % step != 0:
                start = i*L_process + step - ((i*L_process) % step)
            else:
                start = i*L_process

        # > Get sequences for processes
        trunc_seq = seq[start:stop]

        # > Launch task
        job = Process(
            target=worker,
            args=(results, i, k, trunc_seq, full_spec, win, step, fast)
        )
        job.start()
        jobs.append(job)

    # Iterate on tasks --------------------------------------------------------
    for i, job in enumerate(jobs):
        # > Wait the end of the execution
        job.join()

    # > Retrieve results
    for r in results:
        if r is None:
            raise Exception("Error in cumputing, you may reduce number of"
                            "processes")
        distances.extend(r)

    return distances

def worker(results, i, *args, **kwargs):
    results[i] = running_dist(*args, **kwargs)

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

  #######
## TOOLS ######################################################################
  #######

def get_idx(k, seq, win, step):
    return list(range(1, len(seq) - win - k + 3, step))
