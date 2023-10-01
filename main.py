from interpro_json_parser import *
import os
from pymol import cmd, stored
from pyvol import pymol_interface

cmd.reinitialize()

af2_dir = "../../../../../../Hanfeng Modeling PC/AlphaFold2/"
uniprot_id_list = ["A0A0B4J2F2", "O00141", "O00238", "O00311", "O00418", "O00444", "O00506", "O14578", "O14730", "O14733", "O14757", "O14920", "O14936", "O14965", "O14976", "O15021", "O15075", "O15111", "O15146", "O15197"]


interested_ssf = ["SSF56112"]

uniprot_ssf = interpro_json_parser()

files = os.listdir(af2_dir)
pdb_files = [file for file in files if file.endswith('.pdb')]
selected_pdb_files = [pdb_file_name for pdb_file_name in pdb_files if
                      any(uniprot_id in pdb_file_name for uniprot_id in uniprot_id_list)]
print(selected_pdb_files)

# Load all the listed uniprot id into pymol and segment into ssf domains, filter for interested ssf
for file in selected_pdb_files:
    cmd.load(af2_dir + file)
    # pymol_interface.pymol_pocket_cmdline(file.strip(".pdb"))
    uniprotID = file.split("-")[1]
    print(uniprotID)
    for ssf_seq in uniprot_ssf[uniprotID]:
        ssf_name = ssf_seq.pop(0)
        print(ssf_seq)  # this list contains nested lists for domain resID

        selection_all = []
        name_all = [uniprotID]
        for ssf in [[str(element) for element in sublist] for sublist in ssf_seq]:
            [start, end] = ssf
            resID = start + "-" + end
            name_all.append(resID)
            selection = "(" + file.strip(".pdb") + " and resi " + start + "-" + end + ")"
            selection_all.append(selection)
        name_all.insert(1, ssf_name)
        cmd.create('_'.join(name_all), 'or'.join(selection_all))
    cmd.delete(file.strip(".pdb"))
    if interested_ssf:
        for obj in cmd.get_object_list():
            if any(ssf not in obj for ssf in interested_ssf):
                cmd.delete(obj)
    # remove LDDT < 50
    cmd.remove("b<50")

# Call pyvol to calculate the largest pocket in each object


def get_volume(obj):
    with open(obj+".pyvol/"+obj+".rept", encoding="utf-8") as report:
        volume_list = []
        for line in report:
            if line.startswith(obj):
                volume_list.append(float(line.strip().split(",")[1]))
    return volume_list


big_warning_list = []
no_pocket_warning_list = []
for obj in cmd.get_object_list():
    print("[INFO] Calculating for "+obj)
    pymol_interface.pymol_pocket_cmdline(protein=obj, mode="all", min_rad=1.8, max_rad=4, min_volume=500, display_mode="spheres", prefix=obj)
    volumes = get_volume(obj)
    # if no pocket, search min_vol 400?
    if not volumes:
        print("[INFO] Setting min_vol to 400...")
        pymol_interface.pymol_pocket_cmdline(protein=obj, mode="all", min_rad=1.8, max_rad=4, min_volume=400,
                                             display_mode="spheres", prefix=obj)
        volumes = get_volume(obj)
        if not volumes:
            print("[WARNING] No pocket found for "+obj)
            no_pocket_warning_list.append(obj)
    # if volume > 1500, search with max_rad 3.6
    elif volumes[0] > 1500:
        print("[INFO] Setting max_rad to 3.6...")
        big_warning_list.append(obj)
        cmd.delete(obj+"_p*")
        pymol_interface.pymol_pocket_cmdline(protein=obj, mode="all", min_rad=1.8, max_rad=3.6, min_volume=500,
                                             display_mode="spheres", prefix=obj)
        volumes = get_volume(obj)
        # if volume > 1500, search with max_rad 3.2
        if volumes[0] > 1500:
            print("[INFO] Setting max_rad to 3.2...")
            cmd.delete(obj + "_p*")
            pymol_interface.pymol_pocket_cmdline(protein=obj, mode="all", min_rad=1.8, max_rad=3.2, min_volume=500,
                                                 display_mode="spheres", prefix=obj)

print("Too big pocket warning list:" + str(big_warning_list))
print("No pocket warning list:" + str(no_pocket_warning_list))

