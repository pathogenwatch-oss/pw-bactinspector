import json
import subprocess
import sys


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


def get_lineage(taxon_id: str):
    search_str = taxon_id + ','
    with open('data/taxon_names.csv', 'r') as tn:
        for dl in tn.readlines():
            if dl.startswith(search_str):
                names = dl.rstrip().split(',')
                break
    # print(names, file=sys.stderr)

    with open('data/taxon_codes.csv', 'r') as tc:
        for dl in tc.readlines():
            if dl.startswith(search_str):
                codes = dl.rstrip().split(',')
                break
    # print(codes, file=sys.stderr)

    return {
        'names': {'kingdom': names[1], 'genus': names[2], 'species': names[3]},
        'codes': {'kingdom': codes[1], 'genus': codes[2], 'species': codes[3]}
    }


for line in sys.stdin.readlines():
    if line.startswith('file'):
        continue
    if line.rstrip() == '':
        continue
    # print(line, file=sys.stderr)
    data = line.rstrip().split('\t')
    if data[1] == 'No significant matches':
        data = ['', '', '', '', '', '', '1.0', '', '0/1000']
        lineage = {
            'names': {'kingdom': '', 'genus': '', 'species': ''},
            'codes': {'kingdom': '', 'genus': '', 'species': ''}
        }
    else:
        lineage = get_lineage(data[3])
    # if data[3] != '2697049' and data[3] != '90370':
    #     data[2] = data[3]
    # print(data)
    record = {
        'taxId': data[3],
        'speciesId': data[2],
        'speciesName': lineage['names']['species'],
        'genusId': lineage['codes']['genus'],
        "genusName": lineage['names']['genus'],
        "superkingdomId": lineage['codes']['kingdom'],
        "superkingdomName": lineage['names']['kingdom'],
        "referenceId": data[4],
        "mashDistance": data[6],
        "pValue": data[7],
        "matchingHashes": data[8]
    }
    # print(record, file=sys.stderr)
    print(json.dumps(record), file=sys.stdout)
    # break
