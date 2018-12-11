from time import time
from bitstring import Bits
from progressbar import ProgressBar

def get_kmer_freq(k, path):
    start = time()
    seq = str()
    with open(path, 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i == 0:
                continue
            seq += line
    seq = seq.replace('\n', '')
    end_file = time()
    print('length of the sequence:', len(seq))
    d = dict()
    l = len(seq)
    for i in range(l-k):
        kmer = seq[i:i+k]
        if kmer not in d:
            d[kmer] = 1
        else:
            d[kmer] += 1
    end_algo = time()

    print('time file = {:.2e} | time algo = {:.2e}'
            .format(end_file - start, end_algo - end_file))
    return d

def get_kmer_freq_bin(k, path):
    start = time()
    seqs = get_seqs(k, path)
    end_file = time()

    # print('length of the sequence:', len(seq))
    d = dict()
    print('len seqs:', len(seqs))
    bar = ProgressBar()
    for i, seq in bar(enumerate(seqs)):
        start_conv = time()
        N = len(seq)
        seq = convert2bits(seq)
        end_conv = time()
        # print("seq num {} (len={}): {:.2e}s".format(i, N, (end_conv - start_conv)/N))
        l = len(seq)
        for i in range(l-k):
            kmer = seq[i*2:(i+k)*2]
            if kmer not in d:
                d[kmer] = 1
            else:
                d[kmer] += 1
    end_algo = time()

    print('time file = {:.2e} | time algo = {:.2e}'
            .format(end_file - start, end_algo - end_file))
    return d

def load_file(path):
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
    seq = load_file(path)
    seqs = seq.split('N')
    return [s for s in seqs if len(s) >= k]

def main():
    path = './data/GCF_000003645.1_ASM364v1_genomic.fna'
    # for k in (2, 4, int(1e2), int(1e3), int(1e4)):
    #     print('k:', k)
    #     d = get_kmer_freq_bin(k, path)
    get_kmer_freq_bin(int(1e4), path)

# Convertion
str2bit = {
    'A': Bits('0b00'),
    'T': Bits('0b01'),
    'C': Bits('0b10'),
    'G': Bits('0b11')
}
bit2str = {i: k for k,i  in str2bit.items()}

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

if __name__ == "__main__":
    main()
