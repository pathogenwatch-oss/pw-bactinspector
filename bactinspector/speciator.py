import json
import os
import subprocess
import sys

import pandas

from bactinspector import map_to_pw, commands


class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


fasta_path = sys.argv[1]
libraries_path = sys.argv[2]

# print(f'Running {fasta_path} against {libraries_path}\n', file=sys.stderr)
input_dir = os.path.dirname(fasta_path)
fasta_file = os.path.basename(fasta_path)

# Initial mash filter
genus_rep_result = subprocess.run(
    f'mash dist {fasta_path} {libraries_path}/filter.k21s1000.msh -d 0.25 | sort -gk3 | head -n 1 | cut -f 2',
    shell=True, text=True, capture_output=True)

if genus_rep_result.returncode != 0:
    # Search the no genus category
    exit(genus_rep_result.stderr)

genus_rep = genus_rep_result.stdout.replace('./Done/', '').replace('./', '').rstrip()

if genus_rep == '':
    library, threshold = 'NoGenus', 0.05
else:
    print(f'Genus rep: {genus_rep}', file=sys.stderr)
    query_filename = genus_rep.replace('.msh', '.fna.gz')
    taxon_map = pandas.read_parquet(f'{libraries_path}/filter.k21s1000.species.pqt')

    match = taxon_map[taxon_map['filename'] == query_filename].iloc[0]

    if match['superkingdom_name'] != 'Bacteria':
        library, threshold = match['superkingdom_name'], 0.075
    else:
        library, threshold = match['genus_name'].replace(' ', '_'), 0.05

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

json_result = map_to_pw.build_pw_result(results_df)
print(json.dumps(json_result), file=sys.stdout)
