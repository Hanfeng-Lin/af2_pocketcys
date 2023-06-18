import json


def interpro_json_parser():
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
        '''
        print('===============================')
        print('Accession:', accession)
        print('Name:', name)
        print('Source Database:', source_database)
        print('Length:', length)
        print('TaxID:', taxId)
        print('Scientific Name:', scientificName)
        print('Full Name:', fullName)
        print('------------------------------')
        '''

        # Extract entry_subset
        superfamily_list = data['entry_subset']
        ssf_list = []  # [[ssf1,start,end],[ssf2,start,end],...]
        # Iterate over each item in entry_subset
        for ssf in superfamily_list:
            ssf_accession = ssf['accession']
            entry_protein_locations = ssf['entry_protein_locations']

            # Print the extracted entry_subset information
            # print('ssf_Accession:', ssf_accession)

            # Iterate over each entry_protein_locations
            for location in entry_protein_locations:
                fragments = location['fragments']
                model = location['model']
                score = location['score']

                # Print the extracted entry_protein_locations information
                # print('Model:', model)
                # print('Score:', score)

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
                    '''
                    print('Start:', start)
                    print('End:', end)
                    print('DC Status:', dc_status)
                    '''

        uniprot_ssf[accession] = ssf_list

    return uniprot_ssf

