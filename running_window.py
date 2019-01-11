#! /usr/bin/env python3

from pathlib import Path
from time import time

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from progressbar import ProgressBar

from kmerlib.tools import *
from kmerlib.running_window import *
from kmerlib.spectrum import *

# > Seaborn initialization
sns.set()

# Parameters
path = Path('./data/archaea/Aeropyrum_pernix/GCF_000011125.1_ASM1112v1_genomic.fna')
k = 3
win = 3000
step = 100


seq = load_file(path)

# seq = seq[:100000]
# seq = "A"*1013

full_spec = get_kmer_spectrum(k, seq)

start_d = time()
d = running_dist(k, seq, full_spec, win, step)
stop_d = time()
print("Normal: {:.4f}".format(stop_d - start_d))

x = np.array(range(1,1+len(d)*step, step))

plt.figure()
plt.plot(x, d)

# filtered_d = running_average(d, 1)
# plt.figure()
# plt.plot(x, filtered_d, ".-")

start_th = time()
d_th = multithread_running_dist(k, seq, full_spec, win, step, n_thread=2)
stop_th = time()
print("thread: {:.4f}".format(stop_th - start_th))

print("len d", len(d))
print("len d_th", len(d_th))
plt.plot(x, d_th, ".--", label="thread")
plt.legend()

plt.show()

