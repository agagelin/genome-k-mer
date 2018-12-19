from ftplib import FTP
import os
import json
import re
ftp = FTP('ftp.ncbi.nih.gov')

ftp.login()
root = '/genomes/refseq/'
outdir = 'data/'
clades = ['archaea', 'bacteria']

for clade in clades:
    if not os.path.exists(outdir+clade):
        os.makedirs(outdir+clade)
    with open('data/list_orgs/'+clade+'.json') as f:
        ls_orgs_clade = json.load(f)
    for org in ls_orgs_clade:
        print('--------->', org, '<---------')
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
                select_assembly = select_assembly[0] + "/suppressed"
                    
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
                    data = []
                    ftp.dir('-t',data.append)
                    data = [d for d in data if 'GCF_' in d]
                    check_if_assembly = re.search("(GCF_[^\s]+)",data[-1]).group(1)
                if check_if_assembly.__class__ is list:
                    check_if_assembly = check_if_assembly[0]
                ftp.cwd(check_if_assembly)
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
                

                    
                
    
