import re
import sys

import pandas as pd


# header = ',taxid,accession,GCF_accession,GCF_accession_without_version,filename,length,num_contigs,info,refseq_organism_name,refseq_full_organism_name,bacsort_organism_name,curated_organism_name,species_taxid,bioproject,biosample,refseq_category,infraspecific_name,assembly_level,asm_name,submitter,ftp_path,superkingdom_code,genus_code,species_code,superkingdom_name,genus_name,species_name'
# header_fields = header.split(',')


def extract_organism_name(raw: str) -> str:
    subspecies_matcher = re.compile(r'\ssubsp\s\w+$')
    raw = raw.split('/')[0].replace('_', ' ').replace('.fna', '')
    return subspecies_matcher.sub('', raw)


mash_info = {'curated_organism_name': [], 'accession': [], 'length': []}

with open('kleborate_library/mash_info.txt', 'r') as fh:
    for line_number in range(0, 10):
        fh.readline()
    for line in fh.readlines():
        line = line.strip()
        if line == '':
            continue
        data = line.split()
        # Skip comment extension lines for now.
        if data[0] != '1000':
            continue
        mash_info['curated_organism_name'].append(extract_organism_name(data[2]))
        mash_info['accession'].append(data[2])
        mash_info['length'].append(int(data[1]))

mash_df = pd.DataFrame.from_dict(mash_info)
print(mash_df.shape[0], file=sys.stderr)

# Add taxon IDs
taxon_info_df = pd.read_parquet('data/taxon_info.pqt')
print(f'{taxon_info_df.shape[0]} rows in taxon_info df.', file=sys.stderr)
merged = mash_df.merge(taxon_info_df, left_on='curated_organism_name', right_on='species_name', how='inner')
print(f'{merged.shape[0]} rows in merged df.', file=sys.stderr)
merged = merged.drop_duplicates()
print(f'{merged.shape[0]} rows in merged df.', file=sys.stderr)

# Add missing columns
merged['taxid'] = merged['species_code']
merged['species_taxid'] = merged['species_code']
merged['GCF_accession'] = merged['accession']
merged['GCF_accession_without_version'] = merged['accession']
merged['filename'] = merged['accession']
merged['num_contigs'] = 1
merged['info'] = ''
merged['refseq_organism_name'] = merged['curated_organism_name']
merged['refseq_full_organism_name'] = merged['curated_organism_name']
merged['bacsort_organism_name'] = merged['curated_organism_name']
merged = merged.assign(bioproject='', biosample='', refseq_category='', infraspecific_name='', assembly_level='',
                       asm_name='', submitter='', ftp_path='')

# merged = merged.drop(
#     columns=['species_code', 'genus_name', 'genus_code', 'superkingdom_code', 'superkingdom_name']).rename(
#     columns={'species_name': 'organism_name'})

# merged = merged[
#     ['accession', 'GCF_accession', 'GCF_accession_without_version', 'filename', 'length', 'num_contigs', 'info',
#      'refseq_organism_name', 'refseq_full_organism_name', 'bacsort_organism_name', 'curated_organism_name', 'taxid',
#      'species_taxid', 'bioproject', 'biosample', 'refseq_category', 'infraspecific_name', 'assembly_level',
#      'asm_name', 'submitter', 'ftp_path']]

# merged['bioproject'] = ''
# merged['biosample'] = ''
# merged['refseq_category'] = ''

merged.to_parquet('kleborate_library/Kleborate.k21s1000.species.pqt')
merged.to_csv('kleborate_library/Kleborate.k21s1000.species.csv')
