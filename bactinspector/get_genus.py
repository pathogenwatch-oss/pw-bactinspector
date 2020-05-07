import sys

import pandas

query_filename = sys.argv[1].replace('.msh', '.fna.gz').replace('./Done/', '').replace('./', '')
taxon_map = pandas.read_parquet('libraries/filter.k21s1000.species.pqt')
match = taxon_map[taxon_map['filename'] == query_filename].iloc[0]
if match['superkingdom_name'] != 'Bacteria':
    print(match['superkingdom_name'], '0.075', sep='\t', file=sys.stdout)
else:
    print(match['genus_name'].replace(' ', '_'), '0.05', sep='\t', file=sys.stdout)
