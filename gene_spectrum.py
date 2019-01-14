import os
import re
import pandas as pd
from kmerlib.tools import load_seq_file
from kmerlib.spectrum import kmer_spectrum, Spectrum

# def get_coding_sequence (sequence, annotations, concatenated = True):
#     for idx in list(annotations.index):
#         start, stop, name = annotations.loc[idx,["Start", "Stop", "Name"]]

#CDS: Coding Dna Sequence
def compute_cds_spectrum (k, clades = ["archaea", "bacteria"]):
                        #   start_codons = None,
                        #   append_json = True):

#    if not path_csv.endswith(".csv"):
#        raise ValueError("compute_cds_spectrum: file provided was not .csv. Please provide correct file type.")
#    if not path_fna.endswith(".fna"):
#        raise ValueError("compute_cds_spectrum: file provided was not .fna. Please provide correct file type.")
    
    if k>99:
        nk = str(k)
    elif k>9:
        nk = "0" + str(k)
    else:
        nk = "00" + str(k)
    
    path = "data/"
    list_files = list()
    get_file = lambda extention : [f for f  in list_files if f.endswith(extention)]
    
    for clade in clades:
        list_org = os.listdir(path+clade)
        for org in list_org:
            print("#####", org, "#####")
            path_to_org = path+clade+"/"+org+"/"
            list_files = os.listdir(path_to_org)
            csv_file = get_file(".csv")
            fna_file = get_file(".fna")
            # faa_file = get_file(".faa")
            
            whole_spectrum = Spectrum()

            if len(csv_file)>0 and len(fna_file)>0:
                csv_file = csv_file[0]
                fna_file = fna_file[0]

                full_sequence = load_seq_file(path_to_org + fna_file)
                annotations = pd.read_csv(path_to_org + csv_file)
                for idx in list(annotations.index):
                    start, stop = annotations.loc[idx,["Start", "Stop"]]
                    seq = full_sequence[start:stop]
                    seq_spectrum = kmer_spectrum(k, seq, normalize= False)
                    whole_spectrum += seq_spectrum
            
                whole_spectrum.normalize(coef = None)
                json_path = path_to_org + "%s_k%s.json"%(org, nk)
                whole_spectrum.save(json_path)
            else:
                print("Missing fna or csv for", org)




                    
        # annotations = pd.read_csv(path)
        # full_seq = load_file(path_fna)
        

        # x = lambda st, sp: seq[st, sp]
        # pattern_name = re.compile(r".+/(.+)/.+$")
        # pattern_clade = re.compile(r".+/(.+)/.+/.+$")
        # name_org = pattern_name.search(path_csv).group(1)
        # clade = pattern_clade.search(path_csv).group(1)

        # #spectrums = Spectrum()
        # if append_json and not os.path.isfile(json_output):
        #     with open(json_output, "w") as f:
        #         f.close()

        # if start_codons:
        #     print("not implemented")
        
        # else:
        #     for idx in list(annotations.index):
        #         start, stop, name = annotations.loc[idx,["Start", "Stop", "Name"]]

        #         seq = x(start,stop)
            
        #     spectrums[clade][org][name] = get_kmer_spectrum(seq)
        
def main():
    ls_ks = [2, 3, 4, 5, 6, 9]
    ls_ks = [6]
    for k in ls_ks:
        print("\n#########################", k, "#########################")
        compute_cds_spectrum (k = k)

if __name__ == "__main__":
    main()