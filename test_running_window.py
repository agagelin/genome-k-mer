from kmerlib.running_window import *
from kmerlib.spectrum import *
from kmerlib.tools import *

from time import sleep, time

import matplotlib.pyplot as plt

path = "./data/archaea/Aeropyrum_pernix/GCF_000011125.1_ASM1112v1_genomic.fna"

seq = load_file(path)
# seq = seq[:10000]

win = 2000
# Ks = (2, 3, 4, 5, 6, 7, 8, 9, 12, 30, 60, 120)
Ks = (120,)
step = win // 100

# Initialization

# ITERATE ON Ks
for k in Ks:
    t_start = time()
    full_spec = kmer_spectrum(k, seq)
    d = mproc_running_dist(
        k, seq, full_spec, win, step=step, n_process=2)
    t_stop = time()
    sep = " "+FG_GRAY+"|"+END_COLOR+" "
    print(
        FG_GRAY+ "--| " +FG_BLUE+ "DONE"
        +sep+ "K: " +FG_RED+ str(k)
        +sep+ "Win: " +FG_RED+ str(win)
        +sep+ "Step: " +FG_RED+ str(step)
        +sep+ "Time: " +FG_RED+ "{:.3f}".format(t_stop - t_start)
        +END_COLOR+ "sec"
    )

# full_spec = kmer_spectrum(k, seq)
# start = time()
# faster_d = mproc_running_dist(k, seq, full_spec, win, step)
# stop = time()
# print("FASTER", stop - start)

# plt.figure()
# # plt.plot(normal_d)
# plt.plot(faster_d, ".--")
# plt.show()

# print("--> 3 cores")
# ks = (2, 3, 4, 5, 6, 9)
# t = list()
# for k_ in ks:
#     start = time()
#     full_spec = kmer_spectrum(k_, seq)
#     distances = mproc_running_dist(k_, seq, full_spec, win, step, n_process=3)
#     stop = time()
#     print(k_, stop - start)
#     t.append(stop - start)
# plt.plot(ks, t)
# plt.show()

# print("--> 4 cores")
# ks = (2, 3, 4, 5, 6, 9)
# t = list()
# for k_ in ks:
#     start = time()
#     full_spec = kmer_spectrum(k_, seq)
#     distances = mproc_running_dist(k_, seq, full_spec, win, step, n_process=4)
#     stop = time()
#     print(k_, stop - start)
#     t.append(stop - start)
# plt.plot(ks, t)
# plt.show()

# plt.plot(distances)
# plt.show()

# bar = ProgressBar(max_value=10).start()
# for i in range(10):
#     old_val = bar.value
#     bar.update(old_val + 1)
#     sleep(1)
