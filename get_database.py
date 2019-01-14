from ftplib import FTP
import os
import json
import re
import gzip
ftp = FTP('ftp.ncbi.nih.gov')

ftp.login()
root = '/genomes/refseq/'
outdir = 'data/'
clades = ['archaea', 'bacteria']
organism_list = [
    "Halovivax_asiaticus",
    "Acidianus_brierleyi"
]
for clade in clades:
    if not os.path.exists(outdir+clade):
        os.makedirs(outdir+clade)
    with open('data/list_orgs/'+clade+'.json') as f:
        # ls_orgs_clade = json.load(f)

        ls_orgs_clade = organism_list
    for org in ls_orgs_clade:
        print('--------->', org, '<---------', sep = "  ")
        ftp.cwd(root)
        ftp.cwd(clade)
        ls_orgs_ncbi = ftp.nlst()
        current_corresp = []
        i = len(org)
        while len(current_corresp) == 0 and i != 0:
            current_corresp = [ncbi_org for ncbi_org in ls_orgs_ncbi if ncbi_org == org[:i]]
            i -= 1
            
        if i == 0 :
            print('No correspondig organism on ncbi for', org)
        elif len(current_corresp) >1:
            print('More than 1 correspondence for', org)
        else:
            try:
                current_corresp = current_corresp[0]
                if not os.path.exists(outdir+clade+"/"+current_corresp):
                    os.makedirs(outdir+clade+"/"+current_corresp)
                ftp.cwd(current_corresp)
                ls_dir = ftp.nlst()
                select_assembly = [dir for dir in ls_dir if 'latest_assembly' in dir]
                if len(select_assembly) == 0:
                    select_assembly = [dir for dir in ls_dir if 'all_assembly' in dir]
                    if len(select_assembly) > 1:
                        print("More than one dir correspondig to 'all_assembly' for", org, "=> firt selected")
                    select_assembly = select_assembly[0]
                    ftp.cwd(select_assembly)
                    ls_all_assembly = ftp.nlst()
                    if any(i.find("GCF_") >= 0 for i in ls_all_assembly):
                        select_assembly = select_assembly #nothing to do
                    elif any(i == "suppressed" for i in ls_all_assembly):
                        with open(outdir+clade+"/"+current_corresp+"/assembly.suppressed", "w") as f:
                            f.close()
                        print("Only got suppressed assembly for", org)
                        select_assembly = select_assembly + "/suppressed"
                    else:
                        print("No assembly for", org)
                        continue
                    ftp.cwd("..")
                        
                elif len(select_assembly) > 1:
                    print("More than one dir correspondig to 'latest_assembly' for", org, "=> firt selected")
                if select_assembly.__class__ is list:
                    select_assembly = select_assembly[0]
                ftp.cwd(select_assembly)
                ls_assembly = ftp.nlst()
                check_if_assembly = [assembly for assembly in ls_assembly if "GCF_" in assembly]
                if len(check_if_assembly) < 1:
                    print("No assembly for", org)
                else:
                    if len(check_if_assembly) > 1:
                        print("More than one assembly for", org+".", "Taking the last one (the most recent)")
                    entries = list(ftp.mlsd())
                    entries.sort(key = lambda entry: entry[1]['modify'] if "GCF_" in entry[0] else "", reverse = True)
                    right_assembly = entries[0][0]
                    ftp.cwd(right_assembly)
                    fnas = [bin for bin in ftp.nlst() if "genomic.fna.gz" in bin]
                    if len(fnas) < 1:
                        print("No nucleic genome for", org)
                    else:
                        if len(fnas) > 1:
                            fnas = min(fnas, key = len)
                        else:
                            fnas = fnas[0]
                        path_retr = outdir+clade+"/"+current_corresp+"/"+fnas
                        ftp.retrbinary("RETR " + fnas ,open(path_retr, 'wb').write)
                        with gzip.open(path_retr, 'r') as f:
                            file_content = f.read()
                            f.close()
                        with open(path_retr[:-3], "wb") as f:
                            f.write(file_content)
                            f.close()
                        os.remove(path_retr)
                        
                    
                    faa = [bin for bin in ftp.nlst() if "protein.faa.gz" in bin]
                    if len(faa) < 1:
                        print("No protein file for", org)
                    else:
                        if len(faa) > 1:
                            faa = min(faa, key = len)
                        else:
                            faa = faa[0]
                        path_retr = outdir+clade+"/"+current_corresp+"/"+faa
                        ftp.retrbinary("RETR " + faa ,open(path_retr, 'wb').write)
                        with gzip.open(path_retr, 'r') as f:
                            file_content = f.read()
                            f.close()
                        with open(path_retr[:-3], "wb") as f:
                            f.write(file_content)
                            f.close()
                        os.remove(path_retr)
                    

                    gff = [bin for bin in ftp.nlst() if "genomic.gff.gz" in bin]
                    if len(gff) < 1:
                        print("No annotation file for", org)
                    else:
                        if len(gff) > 1:
                            gff = min(gff, key = len)
                        else:
                            gff = gff[0]
                        path_retr = outdir+clade+"/"+current_corresp+"/"+gff
                        ftp.retrbinary("RETR " + gff ,open(path_retr, 'wb').write)
                        with gzip.open(path_retr, 'r') as f:
                            file_content = f.read()
                            f.close()
                        with open(path_retr[:-3], "wb") as f:
                            f.write(file_content)
                            f.close()
                        os.remove(path_retr)
        
            except:
                print("Problem (permission?).")        
                
    
