import json

# The JSON data
#json_data = """[{"metadata":{"accession":"Q13546","name":"Receptor-interacting serine/threonine-protein kinase 1","source_database":"reviewed","length":671,"source_organism":{"taxId":"9606","scientificName":"Homo sapiens","fullName":"Homo sapiens (Human)"}},"entry_subset":[{"accession":"SSF56112","entry_protein_locations":[{"fragments":[{"start":6,"end":286,"dc-status":"CONTINUOUS"}],"model":"0040727","score":7.28e-68}],"protein_length":671,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr011009"},{"accession":"SSF47986","entry_protein_locations":[{"fragments":[{"start":578,"end":666,"dc-status":"CONTINUOUS"}],"model":"0036756","score":2.8e-24}],"protein_length":671,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr011029"}],"taxa":[{"accession":"9606","lineage":["1","131567","2759","33154","33208","6072","33213","33511","7711","89593","7742","7776","117570","117571","8287","1338369","32523","32524","40674","32525","9347","1437010","314146","9443","376913","314293","9526","314295","9604","207598","9605","9606"],"source_database":"uniprot"}]},{"metadata":{"accession":"Q13547","name":"Histone deacetylase 1","source_database":"reviewed","length":482,"source_organism":{"taxId":"9606","scientificName":"Homo sapiens","fullName":"Homo sapiens (Human)"}},"entry_subset":[{"accession":"SSF52768","entry_protein_locations":[{"fragments":[{"start":9,"end":373,"dc-status":"CONTINUOUS"}],"model":"0043894","score":0}],"protein_length":482,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr023696"}],"taxa":[{"accession":"9606","lineage":["1","131567","2759","33154","33208","6072","33213","33511","7711","89593","7742","7776","117570","117571","8287","1338369","32523","32524","40674","32525","9347","1437010","314146","9443","376913","314293","9526","314295","9604","207598","9605","9606"],"source_database":"uniprot"}]}]"""
with open('InterPro_human_superfamily.json', 'r') as file:
    json_data = file.read()

parsed_data = json.loads(json_data)
count=0
# Parse the JSON data
for data in parsed_data:
    count+=1
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

    # Iterate over each item in entry_subset
    for ssf in superfamily_list:
        accession = ssf['accession']
        entry_protein_locations = ssf['entry_protein_locations']

        # Print the extracted entry_subset information
        print('Accession:', accession)

        # Iterate over each entry_protein_locations
        for location in entry_protein_locations:
            fragments = location['fragments']
            model = location['model']
            score = location['score']

            # Print the extracted entry_protein_locations information
            print('Model:', model)
            print('Score:', score)

            # Iterate over each fragment under the same ssf id
            for fragment in fragments:
                start = fragment['start']
                end = fragment['end']
                dc_status = fragment['dc-status']

                # Print the extracted fragment information
                print('Start:', start)
                print('End:', end)
                print('DC Status:', dc_status)

        print()

print(count)