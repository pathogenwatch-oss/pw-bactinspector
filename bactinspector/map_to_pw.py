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


def default_result() -> dict:
    return {
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


def build_pw_result(result_df):
    if result_df.shape[0] == 0:
        return default_result()

    result = result_df.iloc[0]
    species_info = pandas.read_parquet('data/taxon_info.pqt')

    if result['species'] == 'No significant matches':
        return default_result()
    else:
        species_md = species_info.loc[int(result['strain_taxid'])]
        return {
            'taxId': result['strain_taxid'],
            'speciesId': result['species_taxid'],
            'speciesName': species_md['species_name'],
            'genusId': species_md['genus_code'],
            "genusName": species_md['genus_name'],
            "superkingdomId": species_md['superkingdom_code'],
            "superkingdomName": species_md['superkingdom_name'],
            "referenceId": result['top_hit'],
            "mashDistance": float(result['top_hit_distance']),
            "pValue": float(result['top_hit_p_value']),
            "matchingHashes": result['top_hit_shared_hashes']
        }
