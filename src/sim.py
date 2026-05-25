import os
import pandas as pd
import numpy as np
import subprocess

def run_read_analysis(fastq:str, genome_list:str, out_loc:str, threads:int=1):
    
    subprocess.run(['read_analysis.py',
                    'metagenome',
                    '-i', fastq,
                    '-gl', genome_list,
                    '-o', out_loc + '/training/training',
                    '--fastq',
                    '-t', str(threads)], check=True)
    
def run_sim(genome_list:str, abundance_list:str, species_list:str, out_loc:str, perfect:bool=False, threads:int=1):
    
    if perfect:
        subprocess.run(['simulator.py', 'metagenome',
                        '-gl', genome_list,
                        '-a', abundance_list,
                        '-dl', species_list,
                        '-c', out_loc + '/training/training',
                        '-o', out_loc + '/simulated',
                        '--fastq',
                        '--perfect',
                        '-t', str(threads)], check=True)
    else:
        subprocess.run(['simulator.py', 'metagenome',
                        '-gl', genome_list,
                        '-a', abundance_list,
                        '-dl', species_list,
                        '-c', out_loc + '/training/training',
                        '-o', out_loc + '/simulated',
                        '--fastq',
                        '-t', str(threads)], check=True)
        
def run_pbsim_sampling(fastq:str, fasta:str):
    
    subprocess.run(['pbsim',
                    '--strategy', 'wgs',
                    '--method', 'sample',               
                    '--genome', fasta,
                    '--sample', fastq,], check=True)

"""
def run_pbsim_sampling_ext(fastq:str, fasta:str, output:str, depth=None, len_min=None, len_max=None, sample_profile_id=None, 
                           accuracy_min=None, accuracy_max=None, error_ratio=None, bias=None, log_file='pbsim_run_log.txt'):

    prefix = os.path.join(output, 'sd')

    command = [
            'pbsim',
            '--strategy', 'wgs',
            '--method', 'sample',
            '--genome', fasta,
            '--sample', fastq,
            '--prefix', prefix
        ]

    if depth is not None:
        command.extend(['--depth', str(depth)])
    if len_min is not None:
        command.extend(['--length-min', str(len_min)])
    if len_max is not None:
        command.extend(['--length-max', str(len_max)])
    if sample_profile_id is not None:
        command.extend(['--sample-profile-id', str(sample_profile_id)])
    if accuracy_min is not None:
        command.extend(['--accuracy-min', str(accuracy_min)])
    if accuracy_max is not None:
        command.extend(['--accuracy-max', str(accuracy_max)])
    if error_ratio is not None:
        command.extend(['--difference-ratio', str(error_ratio)])
    if bias is not None:
        command.extend(['--hp-del-bias', str(bias)])

    command_str = " ".join(str(x) for x in command)
    print("Running pbsim3 with:", command_str)
    # Run the command and redirect stdout and stderr to a log file only.
    with open(log_file, "a") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
"""

def run_pbsim_sampling_ext(fastq:str, fasta:str, output:str, method:str, depth=None, len_min=None, len_max=None, sample_profile_id=None, 
                           accuracy_min=None, accuracy_max=None, error_ratio=None, bias=None, model=None, len_mean=None, 
                           accuracy_mean=None, log_file='pbsim_run_log.txt'):

    prefix = os.path.join(output, 'sd')

    if model is not None:
        model = "pbsim/" + model + ".model"

    command = [
        'pbsim',
        '--strategy', 'wgs',
        '--method', method,
        '--genome', fasta,
        '--prefix', prefix
    ]

    if method == "sample":
        command.extend(['--sample', fastq])
        """
        command = [
            'pbsim',
            '--strategy', 'wgs',
            '--method', 'sample',
            '--genome', fasta,
            '--sample', fastq,
            '--prefix', prefix
        ]
        """
    elif method == "errhmm":
        command.extend(['--errhmm', model])
    else:
        command.extend(['--qshmm', model])
        """
        command = [
            'pbsim',
            '--strategy', 'wgs',
            '--method', method,
            '--qshmm', model,
            '--genome', fasta,
            '--prefix', prefix
        ]
        """

    if depth is not None:
        command.extend(['--depth', str(depth)])
    if len_min is not None:
        command.extend(['--length-min', str(len_min)])
    if len_max is not None:
        command.extend(['--length-max', str(len_max)])
    if sample_profile_id is not None:
        command.extend(['--sample-profile-id', str(sample_profile_id)])
    if accuracy_min is not None:
        command.extend(['--accuracy-min', str(accuracy_min)])
    if accuracy_max is not None:
        command.extend(['--accuracy-max', str(accuracy_max)])
    if error_ratio is not None:
        command.extend(['--difference-ratio', str(error_ratio)])
    if bias is not None:
        command.extend(['--hp-del-bias', str(bias)])
    if len_mean is not None:
        command.extend(['--length-mean', str(len_mean)])
    if accuracy_mean is not None:
        command.extend(['--accuracy-mean', str(accuracy_mean)])

    command_str = " ".join(str(x) for x in command)
    print("Running pbsim3 with:", command_str)
    # Run the command and redirect stdout and stderr to a log file only.
    with open(log_file, "a") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
