import json

# The JSON data
#json_data = """[{"metadata":{"accession":"Q13546","name":"Receptor-interacting serine/threonine-protein kinase 1","source_database":"reviewed","length":671,"source_organism":{"taxId":"9606","scientificName":"Homo sapiens","fullName":"Homo sapiens (Human)"}},"entry_subset":[{"accession":"SSF56112","entry_protein_locations":[{"fragments":[{"start":6,"end":286,"dc-status":"CONTINUOUS"}],"model":"0040727","score":7.28e-68}],"protein_length":671,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr011009"},{"accession":"SSF47986","entry_protein_locations":[{"fragments":[{"start":578,"end":666,"dc-status":"CONTINUOUS"}],"model":"0036756","score":2.8e-24}],"protein_length":671,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr011029"}],"taxa":[{"accession":"9606","lineage":["1","131567","2759","33154","33208","6072","33213","33511","7711","89593","7742","7776","117570","117571","8287","1338369","32523","32524","40674","32525","9347","1437010","314146","9443","376913","314293","9526","314295","9604","207598","9605","9606"],"source_database":"uniprot"}]},{"metadata":{"accession":"Q13547","name":"Histone deacetylase 1","source_database":"reviewed","length":482,"source_organism":{"taxId":"9606","scientificName":"Homo sapiens","fullName":"Homo sapiens (Human)"}},"entry_subset":[{"accession":"SSF52768","entry_protein_locations":[{"fragments":[{"start":9,"end":373,"dc-status":"CONTINUOUS"}],"model":"0043894","score":0}],"protein_length":482,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr023696"}],"taxa":[{"accession":"9606","lineage":["1","131567","2759","33154","33208","6072","33213","33511","7711","89593","7742","7776","117570","117571","8287","1338369","32523","32524","40674","32525","9347","1437010","314146","9443","376913","314293","9526","314295","9604","207598","9605","9606"],"source_database":"uniprot"}]},{"metadata":{"accession":"A2PYH4","name":"Probable ATP-dependent DNA helicase HFM1","source_database":"reviewed","length":1435,"source_organism":{"taxId":"9606","scientificName":"Homo sapiens","fullName":"Homo sapiens (Human)"}},"entry_subset":[{"accession":"SSF158702","entry_protein_locations":[{"fragments":[{"start":805,"end":986,"dc-status":"CONTINUOUS"}],"model":"0054571","score":1.94e-31}],"protein_length":1435,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":null},{"accession":"SSF46785","entry_protein_locations":[{"fragments":[{"start":688,"end":789,"dc-status":"CONTINUOUS"}],"model":"0054468","score":0.0000169}],"protein_length":1435,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr036390"},{"accession":"SSF52540","entry_protein_locations":[{"fragments":[{"start":293,"end":463,"dc-status":"C_TERMINAL_DISC"},{"start":495,"end":549,"dc-status":"NC_TERMINAL_DISC"},{"start":579,"end":671,"dc-status":"N_TERMINAL_DISC"}],"model":"0052004","score":2.52e-43}],"protein_length":1435,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr027417"}],"taxa":[{"accession":"9606","lineage":["1","131567","2759","33154","33208","6072","33213","33511","7711","89593","7742","7776","117570","117571","8287","1338369","32523","32524","40674","32525","9347","1437010","314146","9443","376913","314293","9526","314295","9604","207598","9605","9606"],"source_database":"uniprot"}]}]"""
with open('InterPro_human_superfamily.json', 'r') as file:
   json_data = file.read()

parsed_data = json.loads(json_data)
count = 0
uniprot_ssf = {}
# Parse the JSON data
for data in parsed_data:
    count += 1
    # Extract metadata
    metadata = data['metadata']
    accession = metadata['accession']
    name = metadata['name']
    source_database = metadata['source_database']
    length = metadata['length']

    # Extract source organism
    source_organism = metadata['source_organism']
    taxId = source_organism['taxId']
    scientificName = source_organism['scientificName']
    fullName = source_organism['fullName']

    # Print the extracted metadata information
    print('===============================')
    print('Accession:', accession)
    print('Name:', name)
    print('Source Database:', source_database)
    print('Length:', length)
    print('TaxID:', taxId)
    print('Scientific Name:', scientificName)
    print('Full Name:', fullName)
    print('------------------------------')

    # Extract entry_subset
    superfamily_list = data['entry_subset']
    ssf_list = []  # [[ssf1,start,end],[ssf2,start,end],...]
    # Iterate over each item in entry_subset
    for ssf in superfamily_list:
        ssf_accession = ssf['accession']
        entry_protein_locations = ssf['entry_protein_locations']

        # Print the extracted entry_subset information
        print('ssf_Accession:', ssf_accession)

        # Iterate over each entry_protein_locations
        for location in entry_protein_locations:
            fragments = location['fragments']
            model = location['model']
            score = location['score']

            # Print the extracted entry_protein_locations information
            print('Model:', model)
            print('Score:', score)

            # Iterate over each fragment under the same ssf id, this means a protein with multiple same domain
            start_end = []
            for fragment in fragments:
                start = fragment['start']
                end = fragment['end']
                dc_status = fragment['dc-status']

                if dc_status == 'CONTINUOUS':
                    ssf_list.append([ssf_accession, [start, end]])
                elif dc_status != 'N_TERMINAL_DISC':
                    start_end.append([start, end])
                else:
                    start_end.append([start, end])
                    ssf_to_append = [ssf_accession]
                    for item in start_end:
                        ssf_to_append.append(item)
                    ssf_list.append(ssf_to_append)





                # Print the extracted fragment information
                print('Start:', start)
                print('End:', end)
                print('DC Status:', dc_status)

    uniprot_ssf[accession] = ssf_list

print(count)
print(uniprot_ssf)
