from interpro_json_parser import *
import os
from pymol import cmd, stored
from pyvol import pymol_interface

uniprot_ssf = interpro_json_parser()

files = os.listdir('./')
pdb_files = [file for file in files if file.endswith('.pdb')]

for file in pdb_files:
    cmd.load(file)
    # pymol_interface.pymol_pocket_cmdline(file.strip(".pdb"))
    uniprotID = file.split("-")[1]
    print(uniprotID)
    for ssf_seq in uniprot_ssf[uniprotID]:
        ssf_seq.pop(0)
        print(ssf_seq)  # this list contains nested lists for domain resID

        selection_all = []
        name_all = [uniprotID]
        for ssf in [[str(element) for element in sublist] for sublist in ssf_seq]:
            [start, end] = ssf
            resID = start+"-"+end
            name_all.append(resID)
            selection = "("+file.strip(".pdb")+" and resi "+start+"-"+end+")"
            selection_all.append(selection)
        cmd.create('_'.join(name_all), 'or'.join(selection_all))
