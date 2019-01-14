import re
import os
import pandas as pd

def annotation_to_csv(path, include_pseudogenes = True):
    """
    Load fna/file.

    Parameter
    ---------
    path: string
        Path to the file

    Output
    ------
    something
    """

    file_content = {"Organism":list(),
                    "Clade":list(),
                    "ID":list(),
                    "Type":list(),
                    "Start":list(),
                    "Stop":list(),
                    "Parent":list(),
                    "Gene ID":list(),
                    "Name":list(),
                    "Gbkey":list(),
                    "Inference":list(),
                    "Product":list(),
                    "Protein ID":list(),
                    "Transl table":list()}

    
    if not path.endswith(".gff"):
        raise ValueError("annotation_to_csv: file provided was not .gff. Please provide correct file type.")
    
    pattern_name = re.compile(r".+/(.+)/.+$")
    pattern_clade = re.compile(r".+/(.+)/.+/.+$")
    name_org = pattern_name.search(path).group(1)
    clade = pattern_clade.search(path).group(1)
    if include_pseudogenes:
        pattern_positions = re.compile(r"^[^#]+\s\w+\s(\w*gene)\s(\d+)\s(\d+).*")
    else:
        pattern_positions = re.compile(r"^[^#^\s]+\s\w+\s(gene)\s(\d+)\s(\d+).*")
    
    pattern_ID = re.compile(r"ID=([^;]+)")
    pattern_parent = re.compile(r"Parent=([^;]+)")
    pattern_geneid = re.compile(r"Dbxref=\w+:([\w\.]+)")
    pattern_name = re.compile(r"Name=([^;]+)")
    pattern_gbkey = re.compile(r"gbkey=([^;]+)")
    pattern_inference = re.compile(r"inference=([^;]+)")
    pattern_product = re.compile(r"product=([^;]+)")
    pattern_proteinid = re.compile(r"protein_id=([^;]+)")
    pattern_table = re.compile(r"transl_table=([^;]+)")
    patterns = [pattern_ID,
                pattern_parent,
                pattern_geneid,
                pattern_name,
                pattern_gbkey,
                pattern_inference,
                pattern_product,
                pattern_proteinid,
                pattern_table]
    keys = ["ID",
            "Parent",
            "Gene ID",
            "Name",
            "Gbkey",
            "Inference",
            "Product",
            "Protein ID",
            "Transl table"]

    detected = False
    with open(path, 'r') as f:
        for line in f:
            if not detected:
                gene_positions = pattern_positions.search(line)
                if gene_positions:
                    detected = True
                    
            else:
                file_content["Organism"].append(name_org)
                file_content["Clade"].append(clade)
                file_content["Type"].append(gene_positions.group(1))
                file_content["Start"].append(int(gene_positions.group(2)) - 1)
                file_content["Stop"].append(int(gene_positions.group(3)) - 1)
                for k, pattern in zip(keys, patterns):
                    result = pattern.search(line)
                    if result:
                        file_content[k].append(result.group(1))
                    else:
                        file_content[k].append(float("nan"))

                detected = False
    
    csv_path = path[:path.rfind("/")]+"/annotation.csv"
    df = pd.DataFrame(file_content)
    df.dropna(axis = 0, how = 'all', inplace=True)
    df.to_csv(csv_path, index = False)




def get_annotation_full_database(include_pseudogenes,
                                 clades = ["archaea", "bacteria"],
                                 skip_suppressed = True):
    path = "data/"
    path = path+"/" if not path.endswith("/") else path
    paths = [path+clade+"/" for clade in clades]
    paths = list()
    for clade in clades:
        list_org = os.listdir(path+clade)
        for org in list_org:
            path_to_gff = path+clade+"/"+org
            list_files = os.listdir(path_to_gff)
            if not ("assembly.suppressed" in list_files and skip_suppressed):
                gff = [f for f in list_files if f.endswith(".gff")]
                if len(gff)>0:
                    gff = gff[0]
                    paths.append(path+clade+"/"+org+"/"+gff)
    
    for path in paths:
        annotation_to_csv(path, include_pseudogenes=include_pseudogenes)


def main():
    get_annotation_full_database(include_pseudogenes = True,
                                 skip_suppressed = False)

if __name__ == "__main__":
    main()