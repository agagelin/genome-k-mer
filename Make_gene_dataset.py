from kmerlib.spectrum import Spectrum, kmer_spectrum
from kmerlib.tools import load_seq_file
from random import sample, seed
from pathlib import Path
import pandas as pd
import json
from itertools import product

def sample_genes_json(genes_per_org,
                n_orgs,
                outfile_name,
                ks = [2,3,4,5,6,9],
                clades = ["archaea", "bacteria"]
                ):
    
    
    result = dict()
    if n_orgs:
        if n_orgs%2 == 0:
            n_archaea = n_orgs/2
            n_bacteria = n_orgs/2
        else:
            n_orgs -= 1
            n_archaea = n_orgs/2
            n_bacteria = n_orgs/2 + 1
        
        paths = list()
        for n, clade in zip([n_archaea, n_bacteria], clades):
            paths.extend(sample(sorted(Path('data').glob('%s/*/'%clade)), n))
        
    else:
        paths = list()
        for clade in clades:
            paths.extend(sorted(Path('data').glob('%s/*/'%clade)))

    for path in paths:
        to_annot = path / 'annotation.csv'
        if to_annot.is_file():
            annotations = pd.read_csv(str(to_annot))
            annotations.dropna(axis = 0, how = 'any', subset = ['Name'], inplace = True)
            to_fna = sorted(path.glob('*.fna'))

            if len(to_fna) > 0 :
                to_fna = to_fna[0]
                seq = load_seq_file(str(to_fna))

                org_name = path.parts[2]
                clade = path.parts[1]
                genes = list(annotations.index)
                if len(genes) >= genes_per_org:
                    kept_genes = sample(genes, genes_per_org)
                    for g in kept_genes:
                        start, stop, name = annotations.loc[g,["Start", "Stop", "Name"]]
                        result[name] = {"Clade":clade,
                                        "Organism":org_name,
                                        "Sequence":seq[start:stop],
                                        "Spectrums":dict()
                                        }
                        for k in ks:
                            result[name]["Spectrums"][k] = kmer_spectrum(k,
                                                                result[name]["Sequence"]).kmer_freqs
                else:
                    print(org_name.replace("_", " "),": too much Nan in gene names")
    outpath = Path("data") / "genes_dataset" / outfile_name
    with open(str(outpath), "w") as file:
        json.dump(result, file, indent=4)
    file.close()

def sample_genes_csv(genes_per_org,
                     n_orgs,
                     set_seed,
                     outfile_name,
                     ks = [2,3,4,5],
                     clades = ["archaea", "bacteria"]
                     ):
    
    
    result = {"Organism" : list(),
              "Clade" : list()
              }
    kmers_dict = dict()
    for k in ks:
        kmers_dict[k] = dict()
        kmers = list(product("ATGC", repeat = k))
        for kmer in kmers:
            kmer_as_string = ''.join(kmer)
            kmers_dict[k][kmer_as_string] = list()
    if n_orgs:
        if n_orgs%2 == 0:
            n_archaea = int(n_orgs/2)
            n_bacteria = int(n_orgs/2)
        else:
            n_orgs -= 1
            n_archaea = int(n_orgs/2)
            n_bacteria = int(n_orgs/2 + 1)

        paths = list()
        for n, clade in zip([n_archaea, n_bacteria], clades):
            if seed:
                seed(set_seed)
            paths.extend(sample(sorted(Path('data').glob('%s/*/'%clade)), n))
        
    else:
        paths = list()
        for clade in clades:
            if seed:
                seed(set_seed)
            paths.extend(sorted(Path('data').glob('%s/*/'%clade)))

    for path in paths:
        to_annot = path / 'annotation.csv'
        if to_annot.is_file():
            annotations = pd.read_csv(str(to_annot))
            annotations.dropna(axis = 0, how = 'any', subset = ['Name'], inplace = True)
            to_fna = sorted(path.glob('*.fna'))

            if len(to_fna) > 0 :
                to_fna = to_fna[0]
                seq = load_seq_file(str(to_fna))

                org_name = path.parts[2]
                clade = path.parts[1]
                print(org_name.replace("_", " "))
                genes = list(annotations.index)
                if len(genes) >= genes_per_org:
                    kept_genes = sample(genes, genes_per_org)
                    for g in kept_genes:
                        start, stop, name = annotations.loc[g,["Start", "Stop", "Name"]]
                        gene_seq = seq[start:stop]
                        result["Clade"].append(clade)
                        result["Organism"].append(org_name)
                        for k in ks:
                            spectrum_dict =  kmer_spectrum(k, gene_seq,
                                                           normalize=True).kmer_freqs
                            for key in kmers_dict[k].keys():
                                if key in spectrum_dict:
                                    kmers_dict[k][key].append(spectrum_dict[key])
                                else:
                                    kmers_dict[k][key].append(0.)

                else:
                    print(org_name.replace("_", " "),": too much Nan in gene names")
    outpath = Path("data") / "genes_dataset" / outfile_name
    # with open(str(outpath), "w") as file:
    #     json.dump(result, file, indent=4)
    # file.close()
    for k in ks:
        result.update(kmers_dict[k])
    df = pd.DataFrame(result)
    df.to_csv(str(outpath), index = False)
    return df
    

def main():
    # sample_genes(genes_per_org = 5,
    #              n_orgs = None,
    #              outfile_name = "5_genes_all_orgs_extended.json")
    sample_genes_csv(genes_per_org = 500,
                     n_orgs = 10,
                     set_seed = 999,
                     outfile_name = "10orgs_500genes_02.csv")

if __name__ == "__main__":
    main()
