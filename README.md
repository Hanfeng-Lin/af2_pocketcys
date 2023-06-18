# af2_pocketcys

something to consider:

too long protein has multiple AF2 files F1, F2... each one with 1400 aa



a domain can be non-continuous

```
{"metadata":{"accession":"Q14156","name":"Protein EFR3 homolog A","source_database":"reviewed","length":821,"source_organism":{"taxId":"9606","scientificName":"Homo sapiens","fullName":"Homo sapiens (Human)"}},"entry_subset":[{"accession":"SSF48371","entry_protein_locations":[{"fragments":[{"start":75,"end":346,"dc-status":"C_TERMINAL_DISC"},{"start":373,"end":493,"dc-status":"N_TERMINAL_DISC"}],"model":"0052026","score":1.55e-12}],"protein_length":821,"source_database":"ssf","entry_type":"homologous_superfamily","entry_integrated":"ipr016024"}],"taxa":[{"accession":"9606","lineage":["1","131567","2759","33154","33208","6072","33213","33511","7711","89593","7742","7776","117570","117571","8287","1338369","32523","32524","40674","32525","9347","1437010","314146","9443","376913","314293","9526","314295","9604","207598","9605","9606"],"source_database":"uniprot"}]}
```



some protein in ssf database are not in human AF proteome (e.g. Q9YNA8). May need call `request` for these.