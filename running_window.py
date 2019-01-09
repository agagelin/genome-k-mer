from pathlib import Path
from kmer_spectrum import get_kmer_spectrum, get_kmer_spectrum_win, load_file
from spectrum_tools import dist, running_dist, multithread_running_dist
from spectrum_tools import running_average
import matplotlib.pyplot as plt
import seaborn as sns
from progressbar import ProgressBar
sns.set()

# Parameters
path = Path('./data/archaea/Aeropyrum_pernix/GCF_000011125.1_ASM1112v1_genomic.fna')
k = 4
win = 3000
step = 100
coef = 1


seq = load_file(path)

# seq = "A"*100000

L = len(seq)
seq = seq[:L//coef]
full_spec = get_kmer_spectrum(k, seq)
d = running_dist(k, seq, full_spec, win, step)

plt.figure()
plt.plot(d)

filtered_d = running_average(d, 100)
plt.figure()
plt.plot(filtered_d)

# d_th = multithread_running_dist(k, seq, full_spec, win, step, n_thread=2)
# print(len(d))
# print(len(d_th))
# plt.figure()
# plt.plot(d_th)

plt.show()

