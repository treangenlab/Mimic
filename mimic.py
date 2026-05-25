''' 
MIMIC Software

@Authors:
Ryan Doughty, Eddy Huang, Yin Min Thant, Iva Kotaskova, Kasambula Arthur Shem, Shwetha Kumar, Mike Nute, Todd Treangen, Mike Nute

@Version 0.2


Mimic creates simulated metagenomes based off of real existing metagenomes.
'''

import argparse
import os
import pathlib
import subprocess
import glob, shutil
import sys
import gzip

from src.tax_identification import run_lemur
from src.sim import run_pbsim_sampling_ext
from src.utils import run_filter, prepare_sample_fasta

__author__ = "Ryan Doughty, Eddy Huang, Yin Min Thant, Iva Kotaskova, Kasambula Arthur Shem, Shwetha Kumar, Mike Nute, Todd Treangen"
__contact__ = "rdd4@rice.edu"

__license__ = "MIT"
__version__ = "0.2"
__email__ = "rdd4@rice.edu"
__status__ = "Development"

def print_info():
    """
    Prints tool information
    """
    print('----------------------------------------')
    print(f'MIMIC Software Version {__version__}')
    print('----------------------------------------')
    print(f'Contact: {__contact__}')
    print(f'Authors: {__author__}\n')
    print(f'Status: {__status__}\n\n')
    
def initialize_working(working:str):
    """Initializes working directory, creates necessary sub directories"""

    if os.path.exists(working):
        raise SystemExit('Initializing Working Directory Failed, Working directory already exists')
    else:
        os.mkdir(working)
        
        lemur = os.path.join(working, 'lemur')
        os.mkdir(lemur)
        
        magnet = os.path.join(working, 'magnet')
        os.mkdir(magnet)

        filter = os.path.join(working, 'filter')
        os.mkdir(filter)
        
        pbsim = os.path.join(working, 'pbsim')
        os.mkdir(pbsim)
        
        simulated_data = os.path.join(working, 'simulated_data')
        os.mkdir(simulated_data)

        print('Initialized Working Directory\n')

def none_if_empty(value):
    """Convert empty string to None."""
    if value is None:
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    return value
    
def run_mimic(args):
    """
    Main pipeline function for the MIMIC
    """
    fastq = args.fastq
    fastq2 = args.fastq2
    output = args.output 
    threads = args.threads
    lemur_db = args.db
    depth = args.depth
    len_min = args.length_min
    len_max = args.length_max
    sample_profile_id = args.sample_profile_id
    accuracy_min = args.accuracy_min
    accuracy_max = args.accuracy_max
    error_ratio = args.difference_ratio
    bias = args.bias
    filter = args.filter
    filter_db = args.filter_db
    method = args.method
    model = args.model
    len_mean = args.length_mean
    accuracy_mean = args.accuracy_mean
    
    pbsim_loc = os.path.join(output, 'pbsim')
    lemur_out = os.path.join(output, 'lemur')
    magnet_out = os.path.join(output, 'magnet')
    filter_out = os.path.join(output, 'filter')

    initialize_working(output)

    if fastq2 is not None:
        combine_out = os.path.join(output, 'combine_out')
        os.mkdir(combine_out)
        combined_fq = os.path.join(combine_out, 'combined')
        combine_cmd = ['pear',
        '-f', fastq,
        '-r', fastq2,
        '-o', combined_fq,
        '-j', str(threads)]
        print("\nRunning pear with:", " ".join(str(x) for x in combine_cmd))
        subprocess.run(combine_cmd, check=True)
        fastq = os.path.join(combine_out, 'combined.assembled.fastq')
    
    ## run lemur
    report_loc = run_lemur(fastq, lemur_db, lemur_out, threads=threads)
    report_loc = os.path.join(lemur_out, 'relative_abundance.tsv')
            
    magnet_cmd = ['python', 'magnet/magnet.py',
                            '-c', report_loc,
                            '-i', fastq,
                            '-o', magnet_out,
                            '-a', '12',
                            '--threads', str(threads)]
    print("\nRunning magnet with:", " ".join(str(x) for x in magnet_cmd))
    subprocess.run(magnet_cmd, check=True)
    
    magnet_report = os.path.join(magnet_out, 'cluster_representative.csv')
    if not os.path.exists(magnet_report):
        raise SystemExit('Magnet failed')
    
    genome_file_loc = os.path.join(magnet_out, 'reference_genomes/merged.fasta')

    if filter:
        genome_file_loc, fastq = run_filter(fastq, magnet_out, filter_out, filter_db, threads)
   
    # Get the current working directory where the files were generated
    working_dir = os.getcwd()
    print("\nPbsim is running from:", working_dir)
    log_file = os.path.join(output, 'pbsim_run_log.txt')

    if(method != "sample"):
        genome_file_loc = prepare_sample_fasta(magnet_out)

    run_pbsim_sampling_ext(fastq, genome_file_loc, pbsim_loc, method, depth, len_min, len_max, sample_profile_id, 
                            accuracy_min, accuracy_max, error_ratio, bias, model, len_mean, accuracy_mean, log_file)
    print("Simulation parameters and log file has been saved to ", output)

    # Decompress all *.fq.gz files in pbsim_loc, concatenate them into one file
    simulated_fastq = os.path.join(output, "simulated.fastq")
    fq_files = glob.glob(os.path.join(pbsim_loc, "*.fq.gz"))

    with open(simulated_fastq, 'wb') as outfile:
        for fq in fq_files:
            print(f"Processing {fq} ...")
            with gzip.open(fq, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

    # Move the concatenated simulated.fastq to the simulated_data directory
    simulated_data_loc = os.path.join(output, 'simulated_data')
    final_path = os.path.join(simulated_data_loc, "simulated.fastq")
    shutil.move(simulated_fastq, final_path)
    print(f"Concatenated fastq file moved to {final_path}")
    
def parse_args():
    parser = argparse.ArgumentParser(description="Universal Taxonomic Classification Verifier.")
    parser.add_argument("-i", "--fastq", type=pathlib.Path, required=True, help="Path to first fastq file")
    parser.add_argument("-I", "--fastq2", type=pathlib.Path, required=False, help="Path to second fastq file for paired-end reads")
    parser.add_argument("-o", "--output", type=pathlib.Path, required=True, help="Path to the output directory.")
    parser.add_argument("--db", type=str, required=True, help='Lemur database location')
    parser.add_argument('-t', '--threads', type=int, required=False, default=1, help='Number of threads for multithreading (Default: 1)')
    parser.add_argument('--depth', type=str, required=False, help='depth of coverage (default 20.0)')
    parser.add_argument('--length-min', type=str, required=False, help='minimum length (default 100)')
    parser.add_argument('--length-max', type=str, required=False, help='maximum length (default 1000000)')
    parser.add_argument('--sample-profile-id', type=str, required=False, help=' sample (filtered) profile ID')
    parser.add_argument('--accuracy-min', type=str, required=False, help='minimum accuracy for simulation, default 0.75')
    parser.add_argument('--accuracy-max', type=str, required=False, help='minimum accuracy for simulation, default 1.00')
    parser.add_argument('--difference-ratio', type=str, required=False, help='difference (error) ratio (default substitution:insertion:deletion = 6:55:39)')
    parser.add_argument('--bias', type=str, required=False, help='bias intensity of deletion in homopolymer')
    parser.add_argument('--filter', action='store_true', help='Will limit simulated reads to user defined databases at --filter_db')
    parser.add_argument('--filter-db', type=str, required=False, help='Path to the filter database')
    parser.add_argument('--method', type=str, required=True, help='pbsim simulating method - sample, errhmm or qshmm')
    parser.add_argument('--model', type=str, required=False, help='pbsim simulating model')
    parser.add_argument('--length-mean', type=str, required=False, help='mean length (default 9000)')
    parser.add_argument('--accuracy-mean', type=str, required=False, help='mean accuracy (default 0.85)')
    

    args = parser.parse_args() 

    # Process optional parameters so that empty strings become None
    args.depth = none_if_empty(args.depth)
    args.length_min = none_if_empty(args.length_min)
    args.length_max = none_if_empty(args.length_max)
    args.sample_profile_id = none_if_empty(args.sample_profile_id)
    args.accuracy_min = none_if_empty(args.accuracy_min)
    args.accuracy_max = none_if_empty(args.accuracy_max)
    args.difference_ratio = none_if_empty(args.difference_ratio)
    args.bias = none_if_empty(args.bias)
    args.model = none_if_empty(args.model)
    args.length_mean = none_if_empty(args.length_mean)
    args.accuracy_mean = none_if_empty(args.accuracy_mean)


    if(args.method != "sample" and args.model is None):
        raise SystemExit("Model must be provided when method is not 'sample'.")
    
    if args.filter and args.filter_db is None:
        raise SystemExit("Filter database must be provided when --filter is set.")
    
    run_mimic(args) 
    

if __name__=='__main__':
    print_info()
    parse_args()
