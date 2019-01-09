class Spectrum:

    # INITIALIZATION ----------------------------------------------------------
    def __init__(self, kmer_freqs=None):
        if kmer_freqs is None:
            self.kmer_freqs = dict()
        else:
            self.kmer_freqs = kmer_freqs

    # Dict behavior -----------------------------------------------------------
    def keys(self):
        return self.kmer_freqs.keys()

    # OPERATOR OVERLOAD ------------------------------------------------------- 
    def __getitem__(self, key):
        return self.kmer_freqs[key]

    def __setitem__(self, key, items):
        self.kmer_freqs[key] = items

    def __contains__(self, key):
        return key in self.kmer_freqs


