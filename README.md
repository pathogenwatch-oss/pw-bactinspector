# Bactinspector
A package to 

1. determine the most probable species based on sequence in fasta/fastq files using refseq and Mash (https://mash.readthedocs.io/en/latest/index.html)
It will count the species of the top ref seq mash matches and report most frequent.
2. determine the closest reference in refseq to a set of fasta/fastq files

The PyPi package is called BactInspectorMax since my original BactInspector asked for the path to a sketch
of refseq genomes. In May 2019 a new curated mash sketch of complete bacterial refseq genome was 
created and bundled into the package. This required a special request to increase the PyPi package limit,
hence the Max.
The command is still `bactinspector <sub command>` however.

## Data
The data bundled into bactinspector are the complete bacterial refseq genomes from May 2019. The species assignations have been corrected by changing species where the majority does not match the species described within the refseq info found [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt) and the excellent [Bacsort resource](https://github.com/rrwick/Bacsort) from [Ryan Wick](https://twitter.com/rrwick)

## Dependencies
[Mash](https://github.com/marbl/Mash/) (> v2.1)
Installation with conda is recommended

## Installation
pip3 install bactinspectorMax

## Usage
```
usage: bactinspector [-h] [-v]
                     {check_species,closest_match,info,create_species_info} ...

A module to determine the most probable species based on sequence in fasta files using refseq and Mash (https://mash.readthedocs.io/en/latest/index.html)
It will count the species of the top ref seq mash matches and report most frequent.

In order to use the module:
  • Specify an input directory and output directory (default is current directory)
  • Specify either a
    • fasta file pattern with -f (e.g "*.fas") or
    • mash sketch file pattern with -m (e.g "*.msh") if you have already sketched the fasta files
  • By default the top 10 matches will be used. Change this with -n
  • Speed things up by changing the number of parallel processes to match the cores on your computer using -p
  • If mash is not in your PATH specify the directory containing the mash executable with -mp

If you want to update the genomes used, follow the instructions on https://gitlab.com/antunderwood/bactinspector/wikis/Updating-the-genomes-in-BactInspector
and use the create_species_info command to make the required file

positional arguments:
  {check_species,closest_match,create_species_info}
                        The following commands are available. Type
                        bactinspector <COMMAND> -h for more help on a specific
                        commands
    check_species       Check the most frequent matches to a species in refseq
    closest_match       Report the closest matches to a set of sequences
    create_species_info
                        Create species info TSV for locally created mash
                        sketches
    info                Provide information about the data in bactinspector
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print out software version
```

### check_species
Assign a species using matches to refseq. Based on observed intra-species mash distances, a result maybe marked as uncertain if the distance to the best hit is greater than 1.2 x the observed maximum intra-species distance.
```
usage: bactinspector check_species [-h] [-i INPUT_DIR] [-o OUTPUT_DIR]
                                   [-p PARALLEL_PROCESSES]
                                   [-n NUM_BEST_MATCHES] [-d DISTANCE_CUTOFF]
                                   [-v ALLOWED_VARIANCE]
                                   [-vl ALLOWED_VARIANCE_RARER_SPECIES] [-s]
                                   [-l LOCAL_MASH_AND_INFO_FILE_PREFIX]
                                   [-mp MASH_PATH]
                                   (-f FASTA_FILE_PATTERN | -fq FASTQ_FILE_PATTERN | -m MASH_SKETCH_FILE_PATTERN)

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        path to input_directory
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        path to output_directory
  -p PARALLEL_PROCESSES, --parallel_processes PARALLEL_PROCESSES
                        number of processes to run in parallel
  -n NUM_BEST_MATCHES, --num_best_matches NUM_BEST_MATCHES
                        number of best matches to return
  -d DISTANCE_CUTOFF, --distance_cutoff DISTANCE_CUTOFF
                        mash distance cutoff (default 0.05)
  -v ALLOWED_VARIANCE, --allowed_variance ALLOWED_VARIANCE
                        proportion of max_distance allowed over which a result
                        will be marked as uncertain (default 0.1)
  -vl ALLOWED_VARIANCE_RARER_SPECIES, --allowed_variance_rarer_species ALLOWED_VARIANCE_RARER_SPECIES
                        proportion of max_distance allowed over which a result
                        will be marked as uncertain for species which have
                        fewer than 10 representatives in refseq (default 0.5)
  -s, --stdout_summary  output a summary of the result to STDOUT
  -l LOCAL_MASH_AND_INFO_FILE_PREFIX, --local_mash_and_info_file_prefix LOCAL_MASH_AND_INFO_FILE_PREFIX
                        the path prefix to the mash sketch file and
                        corresponding info file
  -mp MASH_PATH, --mash_path MASH_PATH
                        path to the mash executable. If not provided it is
                        assumed mash is in the PATH
  -f FASTA_FILE_PATTERN, --fasta_file_pattern FASTA_FILE_PATTERN
                        pattern to match fasta files e.g "*.fas"
  -fq FASTQ_FILE_PATTERN, --fastq_file_pattern FASTQ_FILE_PATTERN
                        pattern to match fastq files e.g "*.fastq.gz"
  -m MASH_SKETCH_FILE_PATTERN, --mash_sketch_file_pattern MASH_SKETCH_FILE_PATTERN
                        pattern to match mash sketch files e.g "*.msh"
```

### closest_match
Find the closest match of a set of genomes to genomes within refseq. Useful as an objective way of choosing a reference genome when mapping

```
usage: bactinspector closest_match [-h] [-i INPUT_DIR] [-o OUTPUT_DIR]
                                   [-p PARALLEL_PROCESSES] [-r]
                                   [-l LOCAL_MASH_AND_INFO_FILE_PREFIX]
                                   [-mp MASH_PATH]
                                   (-f FASTA_FILE_PATTERN | -fq FASTQ_FILE_PATTERN | -m MASH_SKETCH_FILE_PATTERN)

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        path to input_directory
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        path to output_directory
  -p PARALLEL_PROCESSES, --parallel_processes PARALLEL_PROCESSES
                        number of processes to run in parallel
  -r, --ref_and_rep_only
                        only include reference and representative sequences
  -l LOCAL_MASH_AND_INFO_FILE_PREFIX, --local_mash_and_info_file_prefix LOCAL_MASH_AND_INFO_FILE_PREFIX
                        the path prefix to the mash sketch file and
                        corresponding info file
  -mp MASH_PATH, --mash_path MASH_PATH
                        path to the mash executable. If not provided it is
                        assumed mash is in the PATH
  -f FASTA_FILE_PATTERN, --fasta_file_pattern FASTA_FILE_PATTERN
                        pattern to match fasta files e.g "*.fas"
  -fq FASTQ_FILE_PATTERN, --fastq_file_pattern FASTQ_FILE_PATTERN
                        pattern to match fastq files e.g "*.fastq.gz"
  -m MASH_SKETCH_FILE_PATTERN, --mash_sketch_file_pattern MASH_SKETCH_FILE_PATTERN
                        pattern to match mash sketch files e.g "*.msh"
```
### info
Find out what sequences are present in the mash sketch using a term to search the species name. You can either specify `-s` to search the aggregated species information or `-i` to search the individual refseq records.

```
usage: bactinspector info [-h] -t SEARCH_TERM (-s | -i)

optional arguments:
  -h, --help            show this help message and exit
  -t SEARCH_TERM, --search_term SEARCH_TERM
                        search term to use when searching species within
                        bactinspector
  -s, --summary         search the aggregate data
  -i, --individual_records
                        search the individual refseq records
```

### create_species_info

An updated mash sketch file and corresponding info file can be created using the process described [here](https://gitlab.com/antunderwood/bactinspector/wikis/Updating-the-genomes-in-BactInspector)
This process uses the create_species_info commmand

```
usage: bactinspector create_species_info [-h] -m MASH_INFO_FILE -r
                                         REFSEQ_SUMMARY_FILE -b
                                         BACSORT_SPECIES_FILE -x
                                         BACSORT_EXCLUDED_ASSEMBLIES_FILE

optional arguments:
  -h, --help            show this help message and exit
  -m MASH_INFO_FILE, --mash_info_file MASH_INFO_FILE
                        path to info file created using mash info -t
  -r REFSEQ_SUMMARY_FILE, --refseq_summary_file REFSEQ_SUMMARY_FILE
                        path to refseq assembly summary file downloaded via
                        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteri
                        a/assembly_summary.txt
  -b BACSORT_SPECIES_FILE, --bacsort_species_file BACSORT_SPECIES_FILE
                        path to bacsort_species_definitions.txt
  -x BACSORT_EXCLUDED_ASSEMBLIES_FILE, --bacsort_excluded_assemblies_file BACSORT_EXCLUDED_ASSEMBLIES_FILE
                        path to bacsort_excluded_assemblies.txt
```