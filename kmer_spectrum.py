def get_kmer_spectrum(k, full_seq):
    """
    Get the k-mer frequencies.

    Parameters
    ----------
    k: int
        Number of amino-acids.

    full_seq: string
        Sequence.

    Output
    ------
    Spectrum directory built as:
        {k_mer_str: freq, ...}
    """
    seqs = remove_n(k, full_seq)
    spec_dict = dict()
    for _, seq in seqs:
        l = len(seq)
        for i in range(l-k):
            kmer = seq[i:i+k]
            if kmer not in spec_dict:
                spec_dict[kmer] = 1/(l-k)
            else:
                spec_dict[kmer] += 1/(l-k)

    return spec_dict

def get_kmer_spectrum_win(k, full_seq, win):
    """
    Get the k-mer frequencies.

    Parameters
    ----------
    k: int
        Number of amino-acids.

    full_seq: string
        Sequence.

    Output
    ------
    Spectrum directory built as:
        {k_mer_str: freq, ...}
    """
    L = len(full_seq)
    N = L // win
    seqs = remove_n(k, full_seq)
    spec_dict = dict()
    win_specs = [dict() for i in range(N)]
    for n_seq in seqs:
        seq = n_seq[1]
        ns = n_seq[0]
        l = len(seq)
        for i in range(l-k):
            kmer = seq[i:i+k]
            n = ns[i]
            if kmer not in spec_dict:
                spec_dict[kmer] = 1/(l-k)
            else:
                spec_dict[kmer] += 1/(l-k)
            if n <= len(win_specs):
                # print(True)
                if n - win < 0:
                    for spec in win_specs[:n]:
                        if kmer not in spec:
                            spec[kmer] = 1/(l-k)
                        else:
                            spec[kmer] += 1/(l-k)
                else:
                    for spec in win_specs[n - win:n]:
                        print(True)
                        if kmer not in spec:
                            spec[kmer] = 1/win
                        else:
                            spec[kmer] += 1/win

    return spec_dict, win_specs

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

def remove_n(k, seq):
    """
    Split sequences when 'N' is found.
    Remove sequency fragment smaller than k.

    Parameters
    ----------
    k: int
        Number of amino acid in the k-mer.

    seq: string
        Sequence.

    Output
    ------
    Sequencies as list of strings
    """
    seqs = list()
    seqs.append([list(), ''])
    for i, l in enumerate(seq):
        if l == 'N':
            seqs.append([list(), ''])
            continue
        else:
            seqs[-1][0].append(i)
            seqs[-1][1] = seqs[-1][1] + l
    return [s for s in seqs if len(s[1]) >= k]

def main():
    path = './data/GCF_000003645.1_ASM364v1_genomic.fna'
    for k in (2, 4, int(1e2), int(1e3), int(1e4)):
        print('k:', k)
        d = get_kmer_freq(k, path)
    # get_kmer_freq_bin(int(4), path)

if __name__ == "__main__":
    main()
