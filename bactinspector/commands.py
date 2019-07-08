from multiprocessing import Pool
from bactinspector.utility_functions import get_base_name
from bactinspector.mash_functions import run_mash_sketch, get_best_mash_matches, get_most_frequent_species_match, get_species_match_details
import pandas as pd
import datetime
import os, glob, sys



def sample_and_refseq_species_matches(args, num_best_matches):
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

    #  read in species match table
    refseq_species_match_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'all_complete_bacteria_refseq.k21s1000.species.tsv')
    refseq_species_matches = pd.read_csv(refseq_species_match_file, sep = "\t" )
    # filter only reference and representative genomes if specified as an option
    if 'ref_and_rep_only' in args and args.ref_and_rep_only:
        refseq_species_matches = refseq_species_matches.loc[
            (refseq_species_matches['refseq_category'] == 'representative genome') |
            (refseq_species_matches['refseq_category'] == 'reference genome')
        ]

    # run best match processes in parallel
    all_bacterial_refseq_sketches = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'all_complete_bacteria_refseq.k21s1000.msh')
    sample_matches = pool.starmap(get_best_mash_matches, [(sample_sketch,all_bacterial_refseq_sketches, refseq_species_matches, args.output_dir, mash_path, num_best_matches) for sample_sketch in sketch_files])
    return sample_matches, refseq_species_matches
    

def run_check_species(args):
    # get sample matches and refseq_matches
    all_sample_matches, refseq_species_matches = sample_and_refseq_species_matches(args, args.num_best_matches)
        
    results = {'file': [], 'species': [], 'percentage': [], 'top_hit_distance': [], 'top_hit_p_value': [], 'top_hit_shared_hashes': []}
    for filename, sample_matches in all_sample_matches:
        species,  most_frequent_species_count, count, top_hit_distance, top_hit_p_value, top_hit_shared_hashes = get_most_frequent_species_match(sample_matches, refseq_species_matches, args.distance_cutoff)
        results['file'].append(filename)
        results['species'].append(species)
        if species == 'No significant matches':
            results['percentage'].append('NA')
        else:
            results['percentage'].append(int(most_frequent_species_count/count*100))
        results['top_hit_distance'].append(top_hit_distance)
        results['top_hit_p_value'].append(top_hit_p_value)
        results['top_hit_shared_hashes'].append(top_hit_shared_hashes)
    
    results_df = pd.DataFrame(results, columns = ['file', 'species', 'percentage', 'top_hit_distance', 'top_hit_p_value', 'top_hit_shared_hashes']).sort_values('species', ascending=True)
    results_df = results_df.rename(columns = {'percentage' : '%_of_{0}_best_matches=species'.format(args.num_best_matches) })
    now = datetime.datetime.now()
    outfile = os.path.join(args.output_dir, 'species_investigation_{0}.tsv'.format(now.strftime("%Y-%m-%d")))
    results_df.to_csv(outfile, sep = "\t", index = False)
    if args.stdout_summary:
        sys.stdout.write('{0}\n'.format(results_df.to_string(header=False, index=False)))
    sys.stderr.write("Results written to {0}\n".format(outfile))

def run_closest_match(args):
    # get sample matches and refseq_matches
    all_sample_matches, refseq_species_matches = sample_and_refseq_species_matches(args, 1)

    best_matches = None
    for filename, sample_matches in all_sample_matches:
        if best_matches is None:
            best_matches = get_species_match_details(sample_matches, refseq_species_matches)
        else:
            best_matches = best_matches.append(get_species_match_details(sample_matches, refseq_species_matches), )
    
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
