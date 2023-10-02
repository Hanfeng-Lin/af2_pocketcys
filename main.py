from interpro_json_parser import *
import os
import pandas as pd
from pymol import cmd, stored
from pyvol import pymol_interface

cmd.reinitialize()

af2_dir = "../../../../../../Hanfeng Modeling PC/AlphaFold2/"
#uniprot_id_list = ["A0A0B4J2F2", "O00141", "O00238", "O00311", "O00418", "O00444", "O00506", "O14578", "O14730", "O14733", "O14757", "O14920", "O14936", "O14965", "O14976", "O15021", "O15075", "O15111", "O15146", "O15197", "O15264", "O15530", "O43187", "O43283", "O43293", "O43318", "O43353", "O43683", "O43781", "O60229", "O60285", "O60307", "O60566", "O60674", "O75116", "O75385", "O75460", "O75582", "O75676", "O75716", "O75914", "O75962", "O76039", "O94768", "O94804", "O94806", "O94921", "O95382", "O95747", "O95819", "O95835", "O96013", "O96017", "P00519", "P00533", "P00540", "P04049", "P04626", "P04629", "P05129", "P05771", "P06213", "P06239", "P06241", "P06493", "P07332", "P07333", "P07947", "P07948", "P07949", "P08069", "P08581", "P08631", "P08922", "P09619", "P09769", "P0C1S8", "P0C264", "P10398", "P10721", "P11309", "P11362", "P11801", "P11802", "P12931", "P14616", "P15056", "P15735", "P16066", "P16234", "P16591", "P17252", "P17612", "P17948", "P19525", "P19784", "P20594", "P20794", "P21127", "P21675"]
uniprot_id_list= ["P21709", "P21802", "P21860", "P22455", "P22607", "P22612", "P22694", "P23443", "P23458", "P24723", "P24941", "P25092", "P25098", "P27037", "P27361", "P27448", "P28482", "P29317", "P29320", "P29322", "P29323", "P29376", "P29597", "P30291", "P30530", "P31152", "P31749", "P31751", "P32298", "P33981", "P34925", "P34947", "P35590", "P35626", "P35916", "P35968", "P36507", "P36888", "P36894", "P36896", "P36897", "P37023", "P37173", "P41240", "P41279", "P41743", "P42679", "P42680", "P42681", "P42684", "P42685", "P43250", "P43403", "P43405", "P45983", "P45984", "P45985", "P46734", "P48729", "P48730", "P49137", "P49336", "P49674", "P49759", "P49760", "P49761", "P49840", "P49841", "P50613", "P50750", "P51451", "P51617", "P51812", "P51813", "P51817", "P51841", "P51955", "P51956", "P51957", "P52333", "P52564", "P53350", "P53355", "P53667", "P53671", "P53778", "P53779", "P54646", "P54753", "P54756", "P54760", "P54762", "P54764", "P57059", "P57078", "P68400", "P78362", "P78368", "P80192", "Q00526"]
#uniprot_id_list= ["Q00532", "Q00534", "Q00535", "Q00536", "Q00537", "Q01973", "Q01974", "Q02156", "Q02750", "Q02763", "Q02779", "Q02846", "Q04759", "Q04771", "Q04912", "Q05397", "Q05513", "Q05655", "Q05823", "Q06187", "Q06418", "Q07002", "Q07912", "Q08345", "Q08881", "Q09013", "Q12851", "Q12852", "Q12866", "Q13043", "Q13131", "Q13153", "Q13163", "Q13164", "Q13177", "Q13188", "Q13233", "Q13237", "Q13308", "Q13418", "Q13464", "Q13470", "Q13523", "Q13546", "Q13554", "Q13555", "Q13557", "Q13627", "Q13705", "Q13873", "Q13882", "Q13976", "Q14004", "Q14012", "Q14164", "Q14289", "Q14680", "Q15131", "Q15139", "Q15208", "Q15303", "Q15349", "Q15375", "Q15418", "Q15569", "Q15746", "Q15759", "Q15772", "Q15831", "Q15835", "Q16288", "Q16512", "Q16513", "Q16539", "Q16566", "Q16584", "Q16620", "Q16644", "Q16659", "Q16671", "Q16816", "Q16832", "Q2M2I8", "Q32MK0", "Q38SD2", "Q496M5", "Q504Y2", "Q56UN5", "Q59H18", "Q5JZY3", "Q5S007", "Q5TCX8", "Q5TCY1", "Q5VST9", "Q5VT25", "Q6DT37", "Q6IQ55", "Q6J9G0", "Q6P0Q8", "Q6P2M8"]
#uniprot_id_list= ["Q6P3W7", "Q6P5Z2", "Q6PHR2", "Q6SA08", "Q6VAB6", "Q6XUX3", "Q6ZMQ8", "Q6ZN16", "Q6ZWH5", "Q76MJ5", "Q7KZI7", "Q7L7X3", "Q7RTN6", "Q7Z2Y5", "Q7Z7A4", "Q86SG6", "Q86TB3", "Q86TW2", "Q86UE8", "Q86UX6", "Q86V86", "Q86Y07", "Q86YV5", "Q86Z02", "Q8IU85", "Q8IV63", "Q8IVH8", "Q8IVT5", "Q8IVW4", "Q8IW41", "Q8IWB6", "Q8IWQ3", "Q8IWU2", "Q8IY84", "Q8IYT8", "Q8IZE3", "Q8IZL9", "Q8N2I9", "Q8N4C8", "Q8N568", "Q8N5S9", "Q8N752", "Q8NB16", "Q8NCB2", "Q8NE63", "Q8NER5", "Q8NEV1", "Q8NEV4", "Q8NFD2", "Q8NG66", "Q8NI60", "Q8TAS1", "Q8TD08", "Q8TD19", "Q8TDC3", "Q8TDR2", "Q8TDX7", "Q8TEA7", "Q8TF76", "Q8WTQ7", "Q8WU08", "Q8WXR4", "Q8WZ42", "Q92519", "Q92630", "Q92772", "Q92918", "Q96BR1", "Q96C45", "Q96D53", "Q96GD4", "Q96GX5", "Q96J92", "Q96KB5", "Q96KG9", "Q96L34", "Q96L96", "Q96NX5", "Q96PF2", "Q96PN8", "Q96PY6", "Q96Q04", "Q96Q40", "Q96QP1", "Q96QS6", "Q96QT4", "Q96RG2", "Q96RR4", "Q96RU7", "Q96RU8", "Q96S38", "Q96S44", "Q96S53", "Q96SB4", "Q99558", "Q99570", "Q99640", "Q99683", "Q99759", "Q99986"]
#uniprot_id_list= ["Q9BQI3", "Q9BRS2", "Q9BUB5", "Q9BVS4", "Q9BWU1", "Q9BX84", "Q9BXA6", "Q9BXA7", "Q9BXM7", "Q9BXU1", "Q9BYP7", "Q9BYT3", "Q9BZL6", "Q9C098", "Q9C0K7", "Q9H093", "Q9H0K1", "Q9H1R3", "Q9H2G2", "Q9H2K8", "Q9H2X6", "Q9H3Y6", "Q9H422", "Q9H4A3", "Q9H4B4", "Q9H5K3", "Q9H792", "Q9HAZ1", "Q9HBH9", "Q9HBY8", "Q9HC98", "Q9HCP0", "Q9NQU5", "Q9NR20", "Q9NRH2", "Q9NRM7", "Q9NRP7", "Q9NSY1", "Q9NWZ3", "Q9NYL2", "Q9NYV4", "Q9NYY3", "Q9NZJ5", "Q9P0L2", "Q9P1W9", "Q9P286", "Q9P289", "Q9P2K8", "Q9UBE8", "Q9UBS0", "Q9UEE5", "Q9UEW8", "Q9UF33", "Q9UHD2", "Q9UHY1", "Q9UIK4", "Q9UK32", "Q9UKE5", "Q9UKI8", "Q9UL54", "Q9UM73", "Q9UPE1", "Q9UPZ9", "Q9UQ07", "Q9UQ88", "Q9UQB9", "Q9UQM7", "Q9Y243", "Q9Y2H1", "Q9Y2H9", "Q9Y2K2", "Q9Y2U5", "Q9Y3S1", "Q9Y463", "Q9Y4K4", "Q9Y572", "Q9Y5S2", "Q9Y616", "Q9Y6E0", "Q9Y6M4", "Q9Y6R4", "O43930", "P0C263", "P57058", "Q52WX2", "Q5MAI5", "Q6A1A2", "Q7Z695", "Q86YV6", "Q8N165", "Q8NE28", "Q96LW2", "Q9NSY0", "Q9NY57", "Q9Y6S9", "Q3MIX3", "Q6P3R8"]

interested_ssf = ["SSF56112"]

uniprot_ssf = interpro_json_parser()

files = os.listdir(af2_dir)
pdb_files = [file for file in files if file.endswith('.pdb')]
selected_pdb_files = []

# Parsing uniprot_id to find ssf_resID and choose which af2 file should be imported
for uniprot_id in uniprot_id_list:
    pending_pdb_files = []
    for pdb_file_name in pdb_files:
        if uniprot_id in pdb_file_name:
            pending_pdb_files.append(pdb_file_name)  # contains all af2 models for this uniprot ID
    for ssf_seq in uniprot_ssf[uniprot_id]:
        ssf_name = ssf_seq.pop(0)
        print(uniprot_id, ssf_seq, ssf_name)  # ssf_seq contains nested lists for domain resID
        if ssf_name in interested_ssf:
            flattened_resID = []
            [flattened_resID.extend(seg) for seg in ssf_seq]
            # if resID longer than 1400, import Fx model
            if max(flattened_resID) > 1400:
                af2_frag_id = (max(flattened_resID) - 1400)//200 + 2
            else:
                af2_frag_id = 1
            # Import af2 fragment that contains this ssf
            for pdb in pending_pdb_files:
                if "-F"+str(af2_frag_id) in pdb:
                    print("Loading "+pdb)
                    cmd.load(af2_dir + pdb)
                    # re-indexing resID
                    cmd.alter(pdb.strip(".pdb"), "resv += "+str((af2_frag_id-1)*200))
                    # Create object name
                    selection_all = []
                    name_all = [uniprot_id]
                    for ssf in [[str(element) for element in sublist] for sublist in ssf_seq]:
                        [start, end] = ssf
                        resID = start + "-" + end
                        name_all.append(resID)
                        selection = "(" + pdb.strip(".pdb") + " and resi " + start + "-" + end + ")"
                        selection_all.append(selection)
                    name_all.insert(1, ssf_name)
                    print("Creating " + '_'.join(name_all))
                    cmd.create('_'.join(name_all), 'or'.join(selection_all))

                    cmd.delete(pdb.strip(".pdb"))
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
result_dict = {}
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
            volumes = get_volume(obj)

    # return cys within 5A of pocket
    cmd.select("(resn cys and "+obj+") within 5 of "+obj+"_p*")
    stored.cys = set()  # This list will store three ligands: warhead, E3 ligand, and linker.
    cmd.iterate("sele", "stored.cys.add(int(resi))")
    result_dict[obj] = [volumes, len(stored.cys), stored.cys]

print(result_dict)
df = pd.DataFrame(result_dict).transpose()
df.to_csv('uniprot_cys_in_pocket2.csv', index=True, index_label=["uniprot", "pocket_volume", "cys_No.", "cys_resi"])

print("Too big pocket warning list:" + str(big_warning_list))
print("No pocket warning list:" + str(no_pocket_warning_list))

