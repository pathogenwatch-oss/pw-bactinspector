import json
import os
import sys

import pandas

from bactinspector import commands


class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


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


def build_pw_result(result_df, species_datafile):
    if result_df.shape[0] == 0:
        return default_result()

    result = result_df.iloc[0]
    species_info = pandas.read_parquet(species_datafile)

    if result['species'] == 'No significant matches':
        return default_result()
    else:
        species_md = species_info[species_info['species_name'] == result['species']].iloc[0]
        if species_md['species_code'] != result['species_taxid']:
            # This is a re-written reference
            strain_id = species_md['species_code']
        else:
            strain_id = result['strain_taxid']
        return {
            'taxId': strain_id,
            'speciesId': species_md['species_code'],
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


fasta_path = sys.argv[1]
libraries_path = sys.argv[2]

# print(f'Running {fasta_path} against {libraries_path}\n', file=sys.stderr)
input_dir = os.path.dirname(fasta_path)
fasta_file = os.path.basename(fasta_path)

# Initial filter
filter_result = build_pw_result(commands.run_check_species(
    AttributeDict({'distance_cutoff': 0.15,
                   'fasta_file_pattern': fasta_file,
                   'input_dir': input_dir,
                   'output_dir': '/tmp/',
                   'local_mash_and_info_file_prefix': f'{libraries_path}/filter.k21s1000',
                   'stdout_summary': True,
                   'num_best_matches': 10,
                   'parallel_processes': 1,
                   'mash_path': '',
                   'allowed_variance': 0.1,
                   'allowed_variance_rarer_species': 0.5
                   })), 'data/taxon_info.pqt')

genus = filter_result['genusName'].replace(' ', '_') if filter_result['genusName'] != '' else 'NoGenus'

print(f"Genus: {genus}", file=sys.stderr)

escherichia_genus = {'Salmonella', 'Shigella', 'Escherichia'}

if filter_result['superkingdomName'] != 'Bacteria' and filter_result['superkingdomName'] != '':
    library, threshold = filter_result['superkingdomName'], 0.075
elif genus == 'NoGenus':
    library, threshold = genus, 0.075
elif genus in escherichia_genus:
    library, threshold = 'Escherichia', 0.05
else:
    library, threshold = genus, 0.05

print(f'Library={library}, threshold={threshold}', file=sys.stderr)

library_file = f'{libraries_path}/{library}.k21s1000'

results_df = commands.run_check_species(
    AttributeDict({'distance_cutoff': threshold,
                   'fasta_file_pattern': fasta_file,
                   'input_dir': input_dir,
                   'output_dir': '/tmp/',
                   'local_mash_and_info_file_prefix': library_file,
                   'stdout_summary': True,
                   'num_best_matches': 10,
                   'parallel_processes': 1,
                   'mash_path': '',
                   'allowed_variance': 0.1,
                   'allowed_variance_rarer_species': 0.5
                   }))

json_result = build_pw_result(results_df, 'data/taxon_info.pqt')
print(json.dumps(json_result), file=sys.stdout)
