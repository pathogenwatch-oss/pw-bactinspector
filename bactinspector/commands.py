from multiprocessing import Pool
from bactinspector.utility_functions import get_base_name
from bactinspector.mash_functions import run_mash_sketch, get_best_mash_matches, get_most_frequent_species_match, get_species_match_details
from bactinspector.dataframe_parsing_functions  import create_refseq_species_metrics_df, add_certainty_to_merged_results_df, create_refseq_species_info_df
import pandas as pd
import numpy as np
import datetime
import os, glob, sys



def sample_and_refseq_species_info(args, num_best_matches):
    pool=Pool(processes=args.parallel_processes)

    if args.mash_path:
        mash_path = args.mash_path
    else:
        mash_path = ''
    
    if args.fasta_file_pattern:
        # run sketches in parallel
        fasta_files = glob.glob(os.path.join(args.input_dir, args.fasta_file_pattern))
        if len(fasta_files) == 0:
            sys.exit('No files match the pattern {0} in {1}'.format(args.fasta_file_pattern, args.input_dir))
        sketch_files = pool.starmap(run_mash_sketch, [(fasta_file, 'fasta', args.output_dir, mash_path) for fasta_file in fasta_files])
    elif args.fastq_file_pattern:
        fastq_files = glob.glob(os.path.join(args.input_dir, args.fastq_file_pattern))
        if len(fastq_files) == 0:
            sys.exit('No files match the pattern {0} in {1}'.format(args.fastq_file_pattern, args.input_dir))
        sketch_files = pool.starmap(run_mash_sketch, [(fastq_file, 'fastq', args.output_dir, mash_path) for fastq_file in fastq_files])
    elif args.mash_sketch_file_pattern:
        sketch_files = glob.glob(os.path.join(args.input_dir, args.mash_sketch_file_pattern))
        if len(sketch_files) == 0:
            sys.exit('No files match the pattern {0} in {1}'.format(args.mash_sketch_file_pattern, args.input_dir))


    refseq_species_info = create_refseq_species_info_df(args.local_mash_and_info_file_prefix, 'ref_and_rep_only' in args and args.ref_and_rep_only)

    # run best match processes in parallel
    if args.local_mash_and_info_file_prefix is not None:
        all_bacterial_refseq_sketches = args.local_mash_and_info_file_prefix + '.msh'
    else:
        all_bacterial_refseq_sketches = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'all_complete_bacteria_refseq.k21s1000.msh')
    sample_matches = pool.starmap(get_best_mash_matches, [(sample_sketch,all_bacterial_refseq_sketches, refseq_species_info, args.output_dir, mash_path, num_best_matches) for sample_sketch in sketch_files])
    return sample_matches, refseq_species_info
    

def run_check_species(args):
    # get sample matches and refseq_matches
    all_sample_matches, refseq_species_info = sample_and_refseq_species_info(args, args.num_best_matches)
        
    results = {'file': [], 'species': [], 'percentage': [], 'top_hit_distance': [], 'top_hit_p_value': [], 'top_hit_shared_hashes': []}
    for filename, sample_matches in all_sample_matches:
        species,  most_frequent_species_count, count, top_hit_distance, top_hit_p_value, top_hit_shared_hashes = get_most_frequent_species_match(sample_matches, refseq_species_info, args.distance_cutoff)
        results['file'].append(filename)
        results['species'].append(species)
        if species == 'No significant matches':
            results['percentage'].append(np.NAN)
        else:
            results['percentage'].append(int(most_frequent_species_count/count*100))
        results['top_hit_distance'].append(top_hit_distance)
        results['top_hit_p_value'].append(top_hit_p_value)
        results['top_hit_shared_hashes'].append(top_hit_shared_hashes)
    
    # convert to a dataframe and sort
    results_df = pd.DataFrame(results, columns = ['file', 'species', 'percentage', 'top_hit_distance', 'top_hit_p_value', 'top_hit_shared_hashes']).sort_values('species', ascending=True)
    results_df = results_df.rename(columns = {'percentage' : '%_of_{0}_best_matches=species'.format(args.num_best_matches) })

    refseq_species_metrics_df = create_refseq_species_metrics_df()
    results_df = results_df.merge(refseq_species_metrics_df, on = 'species', how = 'left')

    # get df with final result quality of 'good' or 'uncertain'
    results_df = add_certainty_to_merged_results_df(results_df)

    now = datetime.datetime.now()
    outfile = os.path.join(args.output_dir, 'species_investigation_{0}.tsv'.format(now.strftime("%Y-%m-%d")))
    results_df.to_csv(outfile, sep = "\t", index = False)
    if args.stdout_summary:
        sys.stdout.write('{0}\n'.format(results_df.to_string(header=False, index=False)))
    sys.stderr.write("Results written to {0}\n".format(outfile))

def run_closest_match(args):
    # get sample matches and refseq_matches
    all_sample_matches, refseq_species_info = sample_and_refseq_species_info(args, 1)

    best_matches = None
    for filename, sample_matches in all_sample_matches:
        if best_matches is None:
            best_matches = get_species_match_details(sample_matches, refseq_species_info)
        else:
            best_matches = best_matches.append(get_species_match_details(sample_matches, refseq_species_info), )
    
    best_matches['ftp_path'] = best_matches['ftp_path'].str.cat(best_matches['filename'], sep="/")
    best_matches = best_matches.drop(columns = ['distance', 'p-value', 'shared-hashes'])

    # group and sort by count and write to file
    now = datetime.datetime.now()
    outfile = os.path.join(args.output_dir, 'closest_matches_{0}.tsv'.format(now.strftime("%Y-%m-%d")))
    best_matches. \
        fillna(""). \
        groupby(list(best_matches.columns)). \
        size().reset_index(name='count'). \
        sort_values('count', ascending = False). \
        to_csv(outfile, sep = "\t", index = False)
    sys.stderr.write("Results written to {0}\n".format(outfile))

def run_create_species_info(args):

    info_df = pd.read_csv(args.mash_info_file, sep='\t')
    # extract num of contigs
    info_df['num_contigs'] = info_df['Comment'].str.extract(r'\[(\d+) seqs\]')
    # Those rows without [N seqs] have 1 contig
    info_df['num_contigs'].fillna(1, inplace=True)
    # remove [N seqs] from comment and leading spaces
    info_df['Comment'] = info_df["Comment"].str.strip().str.replace('^\[\d+ seqs\] ', '')
    # extract accession and GCF accession
    info_df['accession'] = info_df['Comment'].str.extract(r'^(\w{2}_[\w\.]+)\s')
    info_df['GCF_accession'] = info_df['ID'].str.extract(r'(GCF_\d+\.\d+)')
    info_df['GCF_accession_without_version'] = info_df['GCF_accession'].str.replace(r'\.\d+$', '')
    info_df['Comment'] = info_df["Comment"].str.replace('^\w{2}_[\w\.]+', '').str.replace('\[\.\.\.\]$', '')

    # rename and reorder columns 
    info_df = info_df.rename(columns = {'Length': 'length', 'ID' : 'filename', 'Comment' : 'info'})
    info_df = info_df[['accession', 'GCF_accession', 'GCF_accession_without_version', 'filename', 'length', 'num_contigs', 'info']]

    # read in refseq information
    refseq_assembly_details = pd.read_csv(args.refseq_summary_file, sep="\t", header=1, low_memory = False)
    refseq_assembly_details['GCF_accession_without_version'] = refseq_assembly_details['# assembly_accession'].str.replace(
        r'\.\d+$',
        ''
    )

    # read in bacsort species
    bacsort_species = pd.read_csv(args.bacsort_species_file, sep='\t', comment='#', names = ['GCF_accession_without_version', 'bacsort_organism_name'])
    # read in bacsort excluded
    bacsort_excluded = pd.read_csv(args.bacsort_excluded_assemblies_file, sep='\t', comment='#', names = ['excluded_GCF_accession'])
    
    # merge the mash info and refseq files
    merged = info_df.merge(
        refseq_assembly_details[[
            'GCF_accession_without_version', 'bioproject', 'biosample', 'refseq_category', 'taxid', 'species_taxid',
            'organism_name', 'infraspecific_name', 'assembly_level', 'asm_name', 'submitter', 'ftp_path'
        ]],
        on = ['GCF_accession_without_version'],
        how = 'left'
    )

    # add excluded for those where there are no matches as the most likely cause
    merged.loc[merged['organism_name'].isnull(),'organism_name'] = 'excluded'
    # rename organism_name to full_organism_name
    merged = merged.rename(columns = {'organism_name': 'full_organism_name'})
    merged['organism_name'] = merged['full_organism_name'].str.split(' ').str[0:2].str.join(' ')

    # reorder columns
    merged = merged[[
        'accession', 'GCF_accession', 'GCF_accession_without_version', 'filename', 'length', 'num_contigs', 'info',
        'organism_name', 'full_organism_name', 'taxid', 'species_taxid', 'bioproject', 'biosample','refseq_category',
        'infraspecific_name', 'assembly_level', 'asm_name', 'submitter', 'ftp_path'
    ]]

    # incorporate bacsort 
    merged_with_bacsort =  merged.merge(bacsort_species, on = ['GCF_accession_without_version'], how = 'left')
    merged_with_bacsort =  merged_with_bacsort.merge(bacsort_excluded, left_on = ['GCF_accession_without_version'], right_on = ['excluded_GCF_accession'], how = 'left')

    # curated organism equals bacsort name if it exists
    merged_with_bacsort['curated_organism_name'] = np.where(
        merged_with_bacsort['bacsort_organism_name'].notna(),
        merged_with_bacsort['bacsort_organism_name'],
        np.nan
    )
    merged_with_bacsort['curated_organism_name'] = np.where(
        (merged_with_bacsort['curated_organism_name'].isna() & merged_with_bacsort['excluded_GCF_accession'].notna()),
        'Excluded',
        merged_with_bacsort['curated_organism_name']
    )
    merged_with_bacsort['curated_organism_name'] = np.where(
        merged_with_bacsort['curated_organism_name'].isna(),
        merged_with_bacsort['organism_name'],
        merged_with_bacsort['curated_organism_name']
    )
  
    # convert Salmonella species to serovars
    merged_with_bacsort['salmonella_serovars'] =  merged_with_bacsort['full_organism_name'].str.extract(r'(Salmonella enterica subsp. enterica serovar\s\w+)')
    merged_with_bacsort['curated_organism_name'] = np.where(
        (merged_with_bacsort['full_organism_name'].str.contains('Salmonella enterica subsp. enterica serovar')) & 
        (merged_with_bacsort['curated_organism_name'] == 'Salmonella enterica'),
        merged_with_bacsort['salmonella_serovars'],
        merged_with_bacsort['curated_organism_name'] 
    )
    merged_with_bacsort = merged_with_bacsort.drop(columns = ['salmonella_serovars'])

    # rename and change order of columns
    merged_with_bacsort = merged_with_bacsort.rename(columns = {'organism_name': 'refseq_organism_name', 'full_organism_name': 'refseq_full_organism_name'})
    merged_with_bacsort = merged_with_bacsort[['accession', 'GCF_accession', 'GCF_accession_without_version', 'filename', 'length', 'num_contigs', 'info', 'refseq_organism_name', 'refseq_full_organism_name', 'bacsort_organism_name', 'curated_organism_name', 'taxid', 'species_taxid', 'bioproject', 'biosample','refseq_category', 'infraspecific_name', 'assembly_level', 'asm_name', 'submitter', 'ftp_path']]

    # write out file
    merged_with_bacsort.to_csv(os.path.join(os.path.dirname(args.mash_info_file), 'all_complete_bacteria_refseq.k21s1000.species.tsv'), sep = "\t", index = False)


