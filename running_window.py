from pathlib import Path
from kmer_spectrum import get_kmer_spectrum, get_kmer_spectrum_win, load_file
from spectrum_tools import dist, running_dist
import matplotlib.pyplot as plt
import seaborn as sns
from progressbar import ProgressBar
sns.set()

# Parameters
path = Path('./data/archaea/Aeropyrum_pernix/GCF_000011125.1_ASM1112v1_genomic.fna')
k = 10
win = 200


seq = load_file(path)
full_spec = get_kmer_spectrum(k, seq)
d = running_dist(k, seq, full_spec, win)
print(d)

# d = list()
# bar = ProgressBar()
# for i in bar(range(len(full_seq)-win)):
#     win_spec = get_kmer_spectrum(k, full_seq[i:i+win])
#     d.append(dist(full_spec, win_spec))

# plt.figure()
# plt.plot(d)

