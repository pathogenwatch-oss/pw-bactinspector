import itertools
import subprocess
import sys

import pandas


# Requires the .species.pqt & taxon_info.pqt
# Have a directory with all mash sketches in.
# Writes to a given output directory
def write_metadata(library_name: str, full_df, output_path: str):
    full_df.drop(
        columns=['superkingdom_code', 'superkingdom_name', 'species_code', 'species_name', 'genus_code',
                 'genus_name']).to_parquet(f'{output_path}/{library_name}.pqt')


def write_library(name: str, strain_df, mash_files_path: str, output_path: str):
    write_metadata(name, strain_df, output_path)
    list_filename = f'{name}.lst'
    library_filename = f'{output_path}/{name}.k21s1000.msh'
    strain_names = extract_mash_names(strain_df)
    with open(list_filename, 'w') as list_fh:
        list_fh.write('\n'.join([f'{mash_files_path}/{strainId}.msh' for strainId in strain_names]))
    subprocess.run(f'mash paste -l {list_filename} > {library_filename}', shell=True)


def extract_mash_names(strain_table) -> set:
    return {name.replace('.fna.gz', '') for name in strain_table['filename'].array}


data_dir = sys.argv[1]
mash_dir = sys.argv[2]
output_dir = sys.argv[3]

taxon_data = pandas.read_parquet('data/taxon_info.pqt').rename_axis('taxid')
all_strains = pandas.read_parquet('data/all_complete_refseq.k21s1000.species.pqt').astype({'taxid': int}).set_index(
    'taxid')

# Split by species
viruses = taxon_data[taxon_data['superkingdom_name'] == 'Viruses']
fungi = taxon_data[taxon_data['superkingdom_name'] == 'Eukaryota']
bacteria = taxon_data[taxon_data['superkingdom_name'] == 'Bacteria']

merge_suffixes = ('_strains', '_taxons')
bacteria_strains = all_strains.merge(bacteria, left_on='taxid', right_on='taxid', how='inner', suffixes=merge_suffixes)
virus_strains = all_strains.merge(viruses, left_on='taxid', right_on='taxid', how='inner', suffixes=merge_suffixes)
fungi_strains = all_strains.merge(fungi, left_on='taxid', right_on='taxid', how='inner', suffixes=merge_suffixes)

# Write the virus & fungi libraries
write_library('viruses', virus_strains, mash_dir, output_dir)
write_library('fungi', fungi_strains, mash_dir, output_dir)

# Get the bacteria genus codes
genus_codes = pandas.unique(bacteria['genus_code'])
genus_reps = {genus_code: list() for genus_code in genus_codes}

for genus_code in genus_codes:
    genus_strains = bacteria_strains[bacteria_strains['genus_code'] == genus_code]
    # Write library for all members
    write_library(genus_code, extract_mash_names(genus_strains), mash_dir, output_dir)
    write_metadata(genus_code, genus_strains, output_dir)

    # Get a list of species and pick a representative from each one
    species_codes = pandas.unique(genus_strains['species_code'])
    for species_code in species_codes:
        genus_reps[genus_codes].append(genus_strains[genus_strains['species_code'] == species_code].sample().iloc[0])

# Convert the genus reps into a "strains" dataframe
genus_reps_strains = pandas.DataFrame(list(itertools.chain.from_iterable(genus_reps.values())))

# Create the merged libraries
merged_strains = genus_reps_strains.append(virus_strains).append(fungi_strains)
write_library('filter', merged_strains, mash_dir, output_dir)


# Create a genus-specific library
# strains = taxon_data[taxon_data['genus_code'] == genus_codes]
# with open('genus_code.lst', 'r') as genus_fh:
#     for strain_id in
#         print(strain_id,file=genus_fh)
