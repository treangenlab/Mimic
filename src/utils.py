import os
import argparse
import pathlib
import pandas as pd
import numpy as np
import subprocess
import pysam

from Bio import SeqIO

def parse_metaphlan_file(mpa_report_loc, file_name):
    """
    Parse a MetaPhlAn output file and save the species-level data to a new file for magnet input.
    """
    # Define input and output file paths
    input_file = os.path.join(mpa_report_loc, file_name)
    output_file = os.path.join(mpa_report_loc, "parsed_output.txt")
    
    # Read the input file and skip the first three lines
    with open(input_file, 'r') as fin:
        lines = fin.readlines()[3:]
    
    # Open the output file for writing
    with open(output_file, 'w') as fout:
        # Write the header with tab delimiters
        fout.write("# species_name\tNCBI_tax_id\trelative_abundance\n")
        
        # Process each remaining line
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            # Split the line into fields
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            clade_field = parts[0]
            taxid_field = parts[1]
            rel_abundance = parts[2]
            
            # Split the clade field by '|' to get each taxonomic rank.
            clade_tokens = clade_field.split('|')
            last_token = clade_tokens[-1]
            
            # Only process the line if the last token indicates a species (i.e. starts with "s__")
            if last_token.startswith("s__"):
                species_name = last_token  # Retain the "s__" prefix
                
                # Split the tax_id field and take the last token for the species-level tax_id.
                taxid_tokens = taxid_field.split('|')
                species_taxid = taxid_tokens[-1]
                
                # Write the extracted data to the output file using tab delimiters.
                fout.write(f"{species_name}\t{species_taxid}\t{rel_abundance}\n")
    
    return output_file

def get_abundance(taxid, data):
    try:
        out = data[data['Target_ID'] == taxid]
        if not out.empty:
            return out['F'].values[0]
        else:
            return 0
    except KeyError:
        return 0 

def make_genome_list(metadata_loc, classification_loc, out, working, number_reads):

    metadata = pd.read_csv(metadata_loc)
    classification_data = pd.read_csv(classification_loc, delimiter='\t')
    print(classification_data)

    # get the list for data identified as present in magnet report
    genome_list = metadata[metadata['Presence/Absence'] == 'Present']

    ## get the abundances from classification report and normalize to sum to 100
    genome_list['Abundance'] = genome_list['Taxonomy ID'].apply(lambda row: get_abundance(row, classification_data))

    total_abundance = genome_list['Abundance'].sum()
    if total_abundance > 0:
        genome_list['Abundance'] = np.round(genome_list['Abundance'] / total_abundance * 100, 2)
    
    ## fix the locations names and drop unnecessary taxonomy column
    genome_list['Assembly Accession ID'] = f'{working}/magnet/reference_genomes/' + genome_list['Assembly Accession ID'].astype(str) + '.fasta'

    ## select final columns
    genome_list = genome_list[['Organism of Assembly', 'Assembly Accession ID', 'Abundance']]
    
    abundances =  genome_list[['Organism of Assembly', 'Abundance']] 
    
    out_loc = os.path.join(out, 'abundances.tsv')
    abundances.to_csv(out_loc, sep='\t', header=['Size', number_reads], index=False)
    
    out_loc = os.path.join(out, 'genome_list1.tsv')
    out_loc2 = os.path.join(out, 'genome_list2.tsv')
    genome_list.to_csv(out_loc, sep='\t', header=False, index=False)
    genome_list[['Organism of Assembly', 'Assembly Accession ID']].to_csv(out_loc2, sep='\t', header=False, index=False)
    
    return genome_list

def generate_species_file_info(genome_list, out):
    
    results = []
    for _, row in genome_list.iterrows():
        species = row['Organism of Assembly']
        filename = row['Assembly Accession ID']  
        
        if os.path.exists(filename):
            with open(filename, 'r') as file:
                first_line = file.readline().strip()[1:]
        else:
            first_line = 'File not found'
        
        results.append([species, first_line, 'circular'])
    
    species_info = pd.DataFrame(results, columns=['Species', 'FirstLine', 'Circular'])
    
    out_loc = os.path.join(out, 'species_info.tsv')
    species_info.to_csv(out_loc, sep='\t', header=False, index=False)
    
    return species_info

def run_minimap2(input_fastq, reference_file, output_dir, threads=1):

    print("Run minimap2 to get mapped only sam file")
    mapped_sam_report = os.path.join(output_dir, "aln.sam")

    command = [
        "minimap2", 
        "-ax", "map-ont", 
        str(reference_file), 
        str(input_fastq),
        "-N", str(50),
        "--sam-hit-only",
        "-t", str(threads)
    ]
    
    print("Running minimap2 with:", " ".join(str(x) for x in command))

    with open(mapped_sam_report, "w") as outfile:
        subprocess.run(command, stdout=outfile, stderr=subprocess.DEVNULL, check=True)
    
    return mapped_sam_report

def filter_input_fastq(input_fastq, minimap2_sam, output_dir):
    print("Final filtering")
    filtered_fastq = os.path.join(output_dir, 'filtered.fastq')
    
    # Use pysam to open the SAM file and iterate over alignments
    mapped_reads = set()
    with pysam.AlignmentFile(minimap2_sam, "r") as samfile:
        # Iterate over all alignment records (skips header automatically)
        for record in samfile.fetch(until_eof=True):
            # Add the read name (query_name) to the set
            mapped_reads.add(record.query_name)
    
    # Filter the FASTQ file using Bio.SeqIO based on mapped read names
    with open(input_fastq, "r") as handle, open(filtered_fastq, "w") as output_handle:
        for read in SeqIO.parse(handle, "fastq"):
            if read.id in mapped_reads:
                SeqIO.write(read, output_handle, "fastq")
    
    return filtered_fastq

def run_filter(fastq:str, magnet_out:str, output_dir:str, threads:int=1):
    print("Filtering reads")
    merged_fasta = os.path.join(magnet_out, 'reference_genomes/merged.fasta')
    minimap2_sam = run_minimap2(fastq, merged_fasta, output_dir, threads=threads)

    #filter the original fastq file
    filtered_fastq = filter_input_fastq(fastq, minimap2_sam, output_dir)
    return filtered_fastq

def prepare_sample_fasta(magnet_out):
    sample_fasta = os.path.join(magnet_out, 'sample.fasta')

    cluster_rep = os.path.join(magnet_out, 'cluster_representative.csv')
    df = pd.read_csv(cluster_rep)
    present_df = df[df["Presence/Absence"] == "Present"]
    accession_ids = present_df["Assembly Accession ID"].tolist()
    
    seq_records = []
    for assembly_id in accession_ids:
        reference_fasta = os.path.join(magnet_out, 'reference_genomes', f'{assembly_id}.fasta')

        with open(reference_fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):   
                seq_records.append(record)
                    
    with open(sample_fasta, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")
        
    return sample_fasta

