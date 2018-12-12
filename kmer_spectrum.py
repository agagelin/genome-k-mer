from time import time
from bitstring import Bits
from progressbar import ProgressBar

def get_kmer_spectrum(k, path):
    """
    Get the k-mer frequencies.

    Parameters
    ----------
    k: int
        Number of amino-acids.

    path: string
        Path  to the faa file.

    Output
    ------
    Spectrum directory built as:
        {k_mer_str: freq, ...}
    """

    seqs = get_seqs(k, path)
    spec_dict = dict()
    for seq in seqs:
        l = len(seq)
        for i in range(l-k):
            kmer = seq[i:i+k]
            if kmer not in spec_dict:
                spec_dict[kmer] = 1/(l-k)
            else:
                spec_dict[kmer] += 1/(l-k)

    return spec_dict

def get_kmer_freq_bin(k, path):
    """
    Get the k-mer frequencies with k-mers converted in binaries.

    Parameters
    ----------
    k: int
        Number of amino-acids.

    path: string
        Path  to the faa file.

    Output
    ------
    Spectrum directory built as:
        {k_mer_str: freq, ...}
    """
    seqs = get_seqs(k, path)
    spec_dict = dict()
    bar = ProgressBar()
    for i, seq in bar(enumerate(seqs)):
        N = len(seq)
        seq = convert2bits(seq)
        l = len(seq)
        for i in range(l-k):
            kmer = seq[i*2:(i+k)*2]
            if kmer not in spec_dict:
                spec_dict[kmer] = 1/(l-k)
            else:
                spec_dict[kmer] += 1/(l-k)

    return spec_dict

def load_file(path):
    """
    Load fna/file.

    Parameter
    ---------
    path: string
        Path to the file

    Output
    ------
    Sequence as a string.
    """
    seq = str()
    with open(path, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i == 0:
                continue
            seq += line
    seq = seq.replace('\n', '')
    print('Length of the genome:', len(seq))
    return seq

def get_seqs(k, path):
    """
    Split sequences when 'N' is found.
    Remove sequency fragment smaller than k.

    Parameters
    ----------
    k: int
        Number of amino acid in the k-mer.

    path: string
        Path to the file

    Output
    ------
    Sequencies as list of strings
    """
    seq = load_file(path)
    seqs = seq.split('N')
    return [s for s in seqs if len(s) >= k]

def convert2bits(string):
    seq = Bits()
    for letter in string:
        seq += str2bit[letter]
    return seq

def convert2str(bits):
    seq = ""
    N = int(len(bits)/2)
    for i in range(N):
        seq += str2bit[bits[i*2:i*2+2]]
    return seq

def main():
    path = './data/GCF_000003645.1_ASM364v1_genomic.fna'
    for k in (2, 4, int(1e2), int(1e3), int(1e4)):
        print('k:', k)
        d = get_kmer_freq(k, path)
    # get_kmer_freq_bin(int(4), path)

# Convertion
str2bit = {
    'A': Bits('0b00'),
    'T': Bits('0b01'),
    'C': Bits('0b10'),
    'G': Bits('0b11')
}
bit2str = {i: k for k,i  in str2bit.items()}

if __name__ == "__main__":
    main()
