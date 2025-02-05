# Import packages.
import gzip
import numpy as np
import pandas as pd
import sys

### sys.argv[1]: path to the passed qc table with the suffix txt.gz ###
### sys.argv[1]: path to the mono- and bi-allelic hla pbd cds vcf with the suffix txt.gz ###
### sys.argv[2]: path to the vcf.gz file ###


# Define a function to filter the original vcf file.
def filter_replaced_vcf(qced_file, hla_pbd_cds_vcf, vcf, buffer_size=100_000):
    """
    Outputs the filtered VCF file to stdout.
    """
    # Load the qc'ed file.
    qced_df = pd.read_csv(qced_file, sep='\t')
    # Extract the qc ids and positions.
    qc_ids = qced_df.qc_id.values
    all_qced_pos = qced_df.pos.values
    # Intialize a set to store the qc'ed postions.
    qced_pos = set(all_qced_pos[qc_ids == 0])
    
    # Load the hla pbd cds vcf.
    hla_pbd_cds_df = pd.read_csv(hla_pbd_cds_vcf, sep='\t')
    # Extract the alternative alleles.
    hla_pbd_cds_alts = hla_pbd_cds_df.ALT.values
    # Intialize allelic masks.
    is_mono = hla_pbd_cds_alts == '.'
    is_bi = np.array([len(alt) == 1 for alt in hla_pbd_cds_alts]) & ~is_mono
    # Filter the hla pbd cds vcf only for bi-allelic sites.
    hla_pbd_cds_vcf_df = hla_pbd_cds_df[is_bi].reset_index(drop=True)
    # Extract the hla pbd cds positions.
    hla_pbd_cds_pos = hla_pbd_cds_vcf_df.POS.values
    
    # Intialize a list to store the qc'ed vcf lines.
    vcf_lines = []
    
    # Read the original vcf file.
    with gzip.open(vcf, 'rt') as vcf_data:
        # Iterate through every line in the original vcf file.
        for line in vcf_data:
            # If the line contains meta information.
            if line.startswith('#'):
                # Append the line to the vcf list.
                vcf_lines.append(line)
            # Else-if the current position is a new or replaced hla pbd cds call.
            elif int(line.split()[1]) in hla_pbd_cds_pos:
                # Append the line to the vcf list.
                vcf_lines.append('\t'.join(hla_pbd_cds_vcf_df[hla_pbd_cds_pos == int(line.split()[1])].values[0].astype(str))+'\n')
            # Else-if the current line passed qc.
            elif int(line.split()[1]) in qced_pos:
                # Append the line to the vcf list.
                vcf_lines.append(line)
            
            # If the number of vcf lines exceeds the buffer size.
            if len(vcf_lines) >= buffer_size:
                # Write the vcf lines to the stdout.
                sys.stdout.writelines(vcf_lines)
                # Clear the written vcf lines.
                vcf_lines.clear()

        # If there are still vcf lines to be written.
        if vcf_lines:
            # Write the remaining vcf lines to stdout.
            sys.stdout.writelines(vcf_lines)
    return

# Filter the replaced vcf file.
filter_replaced_vcf(str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]))