from progressbar import ProgressBar

def get_kmer_spectrum(k, seq):
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
    spec_dict = dict()
    L = len(seq)
    for i in range(L-k):
        kmer = seq[i:i+k]
        if "N" not in kmer:
            if kmer not in spec_dict:
                spec_dict[kmer] = 1/(L-k)
            else:
                spec_dict[kmer] += 1/(L-k)

    return spec_dict

def get_kmer_spectrum_win(k, seq, win):
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
    spec_dict = dict()
    L = len(seq)
    win_specs = [dict() for i in range(L-win)]
    bar = ProgressBar()
    for i in bar(range(L-k)):
        kmer = seq[i:i+k]
        if "N" in kmer:
            continue
        if kmer not in spec_dict:
            spec_dict[kmer] = 1/(L-k)
        else:
            spec_dict[kmer] += 1/(L-k)
        if i <= len(win_specs):
            if i - win < 0:
                for spec in win_specs[:i]:
                    if kmer not in spec:
                        spec[kmer] = 1/(L-k)
                    else:
                        spec[kmer] += 1/(L-k)
            else:
                for spec in win_specs[i - win:i]:
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
    seqs = seq.split('N')
    return [s for s in seqs if len(s[1]) >= k]

def main():
    path = './data/GCF_000003645.1_ASM364v1_genomic.fna'
    for k in (2, 4, int(1e2), int(1e3), int(1e4)):
        print('k:', k)
        seq = load_file(path)
        d = get_kmer_spectrum(k, seq)
        print(d)
    # get_kmer_freq_bin(int(4), path)

if __name__ == "__main__":
    main()
