import json
import argparse
import os
from time import time
from pathlib import Path

import numpy as np

from kmerlib.running_window import *
from kmerlib.spectrum import *
from kmerlib.tools import *
from utils.term_colors import *

# ARGPARSER -------------------------------------------------------------------
parser = argparse.ArgumentParser()
# Input file
parser.add_argument("organism", type=str,
                    help="Organism partial path. Example: 'archaea/Aeropyrum'")
# K
parser.add_argument("-k", type=int, nargs="+",
                    default=(2, 3, 4, 5, 6, 9, 12, 24, 30, 60, 120),
                    help="Number of nucleic acids in the K-mer.")
# Window length
parser.add_argument("-w", "--win", type=int, nargs="+",
                    default=(200, 500, 1000, 2000, 5000, 10000, 50000, 100000),
                    help="Size of the window")
# Output file
parser.add_argument("-o", "--output", type=str,
                    default=None,
                    help="Output file path.")
# Number of process
parser.add_argument("-n", "--n-process", type=int, default=3,
                    help="number of process to run")

script_args = parser.parse_args()

# Load data -------------------------------------------------------------------

path = Path("./data") / script_args.organism
files = os.listdir(str(path))
fna_files = [f for f in files if ".fna" in f]
if len(fna_files) != 1:
    raise Exception("Several .fna files or no .fna file")
path = path / fna_files[0]

type_, name = script_args.organism.split("/")
data = {
    "name": name,
    "type": type_,
    "path": str(path),
    "seq": load_seq_file(str(path))
}

# Settings --------------------------------------------------------------------
Ks = script_args.k
win_lens = script_args.win

# Computation -----------------------------------------------------------------

# Initialization
data["distances"] = dict()
seq = data["seq"]
n_process = script_args.n_process
sep = " "+FG_GRAY+"|"+END_COLOR+" "
print("> " +FG_GREEN+ name +sep+FG_BLUE+ type_ +END_COLOR)

# Iterate on Ks
for k in Ks:

    data["distances"][k] = dict()
    full_spec = kmer_spectrum(k, seq)

    # Iterate on Window length
    for win_len in win_lens:
        step = win_len // 100

        t_start = time()
        data["distances"][k][win_len] = mproc_running_dist(
            k, seq, full_spec, win_len, step=step, n_process=n_process
        )
        t_stop = time()

        print(
            FG_GRAY+ "--| " +FG_BLUE+ "DONE"
            +sep+ "K: " +FG_GREEN+ str(k)
            +sep+ "Win: " +FG_GREEN+ str(win_len)
            +sep+ "Step: " +FG_GREEN+ str(step)
            +sep+ "Time: " +FG_GREEN+ "{:.2f}".format(t_stop - t_start)
            +END_COLOR+ "sec"
        )
    
# Save data
output_path = script_args.output
if output_path is None:
    file_name = "{}_k{}-{}_w{}-{}.json".format(
                 name, Ks[0], Ks[-1], win_lens[0], win_lens[-1]
    )
    output_path = Path("./data/running_window_data") / file_name
else:
    output_path = Path(output_path)

with open(str(output_path), "w") as f:
    json.dump(data, f)
print(FG_GRAY+ "--| " +END_COLOR+ "Written as: " +FG_YELLOW+ str(output_path)
      +END_COLOR)
