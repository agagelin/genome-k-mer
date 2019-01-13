from time import time
import json

import numpy as np

from kmerlib.running_window import *
from kmerlib.spectrum import *
from kmerlib.tools import *
from utils.term_colors import *

# Load data -------------------------------------------------------------------

data = [dict() for i in range(4)]

data[0]["name"] = "Aeropyrum permnix"
data[0]["type"] = "archaea"
data[0]["path"] = "./data/archaea/Aeropyrum_pernix/GCF_000011125.1_ASM1112v1_genomic.fna"

data[1]["name"] = "Methanocaldococcus fervens"
data[1]["type"] = "archaea"
data[1]["path"] = "./data/archaea/Methanocaldococcus_fervens/GCF_000023985.1_ASM2398v1_genomic.fna"

data[2]["name"] = "Acholeplasma laidlawii"
data[2]["type"] = "bacteria"
data[2]["path"] = "./data/bacteria/Acholeplasma_laidlawii/GCF_900476025.1_50465_E02_genomic.fna"

data[3]["name"] = "Campylobacter jejuni"
data[3]["type"] = "bacteria"
data[3]["path"] = "./data/bacteria/Campylobacter_jejuni/GCF_900638365.1_57428_D01_genomic.fna"

for s in data:
    s["seq"] = load_file(s["path"])

# Settings --------------------------------------------------------------------
Ks = (2, 3, 4, 5, 6, 9, 12, 24)
win_lens = (200, 500, 1000, 2000, 5000, 10000, 50000, 100000)

# Computation -----------------------------------------------------------------
# Iterate on species
for s in data:
    # Initialization
    s["distances"] = dict()
    name = s["name"]
    type_ = s["type"]
    seq = s["seq"]
    print("> " +FG_RED+ name +" "+FG_BLUE+ type_ +END_COLOR)
    
    # ITERATE ON Ks
    for k in Ks:
        if k >= 60:
            n_process = 1
        else:
            n_process = 3

        s["distances"][k] = dict()
        full_spec = kmer_spectrum(k, seq)

        # Iterate on Window length
        for win_len in win_lens:
            step = win_len // 100

            t_start = time()
            s["distances"][k][win_len] = mproc_running_dist(
                k, seq, full_spec, win_len, step=step, n_process=n_process
            )
            t_stop = time()

            sep = " "+FG_GRAY+"|"+END_COLOR+" "
            print(
                FG_GRAY+ "--| " +FG_BLUE+ "DONE"
                +sep+ "K: " +FG_RED+ str(k)
                +sep+ "Win: " +FG_RED+ str(win_len)
                +sep+ "Step: " +FG_RED+ str(step)
                +sep+ "Time: " +FG_RED+ "{:.3f}".format(t_stop - t_start)
                +END_COLOR+ "sec"
            )
    
# Save data
with open("data/running_window_data.json", "w") as f:
    json.dump(data, f)
