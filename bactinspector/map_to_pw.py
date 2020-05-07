import json
import sys
import pandas

# {
#     "taxId": "90370",
#     "speciesId": "28901",
#     "speciesName": "Salmonella enterica",
#     "genusId": "590",
#     "genusName": "Salmonella",
#     "superkingdomId": "2",
#     "superkingdomName": "Bacteria",
#     "referenceId": "GCF_001084525.1",
#     "mashDistance": 0.00015683822533911295,
#     "pValue": 0,
#     "matchingHashes": "398/400"
# }

species_info = pandas.read_parquet('data/taxon_info.pqt')
for line in sys.stdin.readlines():
    if line.startswith('file'):
        continue
    if line.rstrip() == '':
        continue
    data = line.rstrip().split('\t')
    if data[1] == 'No significant matches':
        record = {
            'taxId': '',
            'speciesId': '',
            'speciesName': '',
            'genusId': '',
            "genusName": '',
            "superkingdomId": '',
            "superkingdomName": '',
            "referenceId": '',
            "mashDistance": 1.0,
            "pValue": '',
            "matchingHashes": '0/1000'
        }
    else:
        species_md = species_info.loc[int(data[3])]
        record = {
            'taxId': data[3],
            'speciesId': data[2],
            'speciesName': species_md['species_name'],
            'genusId': species_md['genus_code'],
            "genusName": species_md['genus_name'],
            "superkingdomId": species_md['superkingdom_code'],
            "superkingdomName": species_md['superkingdom_name'],
            "referenceId": data[4],
            "mashDistance": float(data[6]),
            "pValue": float(data[7]),
            "matchingHashes": data[8]
        }
    print(json.dumps(record), file=sys.stdout)
    break
