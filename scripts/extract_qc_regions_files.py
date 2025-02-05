# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1]: dataset prefix ###


# Define a function to load the qc table.
def load_qc_table(prefix):
    """
    Returns the vcf qc table as pandas dataframe and dictionary.
    """
    df = pd.read_csv(f'./meta_data/{prefix}_qc_chr6.txt.gz', sep='\t')
    dicc = {col: df[col].values for col in df.columns.values}
    return df, dicc


# Define a function to extract the regions files.
def extract_regions_file(prefix):
    """
    Returns the regions files for the sites that passed and failed qc.
    """
    # Load the qc table.
    org_qc_df, org_qc_info = load_qc_table(f'{prefix.lower()}_original')
    # Intialize the qc'ed mask.
    is_qced = org_qc_info['qc_id'] == 0
    # Export the qc results.
    org_qc_df[is_qced].to_csv(f'./meta_data/{prefix.lower()}_original_passed_qc_chr6.txt.gz', sep='\t', index=False)
    org_qc_df[~is_qced].to_csv(f'./meta_data/{prefix.lower()}_original_failed_qc_chr6.txt.gz', sep='\t', index=False)
    return

# Extract the region files.
extract_regions_file(str(sys.argv[1]))