import pandas as pd
import numpy as np
import os


def create_refseq_species_info_df(local_mash_and_info_file_prefix, ref_and_rep_only = False):
    #  read in species match table
    if local_mash_and_info_file_prefix is not None:
        refseq_species_info_file = local_mash_and_info_file_prefix + '.species.tsv'
    else:
        refseq_species_info_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'all_complete_bacteria_refseq.k21s1000.species.tsv')
    refseq_species_info = pd.read_csv(refseq_species_info_file, sep = "\t" )


    # filter only reference and representative genomes if specified as an option
    if ref_and_rep_only:
        refseq_species_info = refseq_species_info.loc[
            (refseq_species_info['refseq_category'] == 'representative genome') |
            (refseq_species_info['refseq_category'] == 'reference genome')
        ]
    return refseq_species_info
    

def create_refseq_species_metrics_df():
    """
    Create a pandas dataframe based on the per species distance metrics found in the file
    all_complete_bacteria_refseq.k21s1000.species.distance_and_length_metrics.tsv
    """
    # merge with species distances cutoff and length
    refseq_species_metrics_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'all_complete_bacteria_refseq.k21s1000.species.distance_and_length_metrics.tsv')
    refseq_species_metrics_df = pd.read_csv(refseq_species_metrics_file, sep = '\t')
    


    # add adjusted distance 0.1 for max_distance > 0.1
    refseq_species_metrics_df['adjusted_max_distance'] = np.where(refseq_species_metrics_df['max_distance'] > 0.1, 0.1, refseq_species_metrics_df['max_distance'] )

    # where there are few distances: replace max_length and min_length with 'unknown'
    refseq_species_metrics_df['max_length'] = np.where(refseq_species_metrics_df['max_distance'] == 0, np.nan, refseq_species_metrics_df['max_length'])
    refseq_species_metrics_df['min_length'] = np.where(refseq_species_metrics_df['max_distance'] == 0, np.nan, refseq_species_metrics_df['min_length'])

    # fill NaNs with 'Unknown'
    refseq_species_metrics_df = refseq_species_metrics_df.fillna('unknown')
    return refseq_species_metrics_df


def add_certainty_to_merged_results_df(results_df, distance_threshold_extension = 1.1, num_best_matches = 10):
    results_df['result'] = np.where(
        (results_df['max_distance'] == 0 ) |
        (results_df['max_distance'].isna()  & results_df['top_hit_distance'] > 0.05) |
        (results_df['top_hit_distance'] > results_df['adjusted_max_distance'] * distance_threshold_extension) |
        (results_df['%_of_{0}_best_matches=species'.format(num_best_matches)].isna()) |
        (results_df['%_of_{0}_best_matches=species'.format(num_best_matches)] < 60),
        'uncertain', 'good'
    )

    # convert lengths to int
    results_df['max_length'] = pd.to_numeric(results_df['max_length'], errors = 'ignore', downcast='integer')
    results_df['min_length'] = pd.to_numeric(results_df['min_length'], errors = 'ignore', downcast='integer')

    # rename lengths
    results_df = results_df.rename(
        columns = {
            'max_length' : 'maximum_genome_length',
            'min_length' : 'minimmum_genome_length'
        }
    )
    # drop columns no longer required
    results_df = results_df.drop(columns = [
        'max_distance',
        'mean_distance',
        'std_distance',
        'quartile_75_distance',
        'quartile_25_distance',
        'adjusted_max_distance']
        )
    return results_df