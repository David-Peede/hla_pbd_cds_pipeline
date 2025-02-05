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

# Define a function to compute effective sequence lengths.
def compute_eff_seq_lens(prefix):
    """
    Returns the effective sequence lengths for the exons and focal regions as .csv.gz .
    """
    # Load the exons and focal regions.
    exons_df = pd.read_csv(
        './meta_data/ncbi_refseq_select_chr6_exons.txt',
        usecols=[1, 2], names=['start', 'end'], sep='\t'
    )
    regions_df = pd.read_csv(
        './meta_data/chr6_focal_regions.txt',
        names=['region', 'start', 'end'], sep='\t'
    )
    # Compute the length of each exon and region.
    exons_df['length'] = exons_df['end'] - exons_df['start'] + 1
    regions_df['length'] = regions_df['end'] - regions_df['start'] + 1

    # Load the failed qc information.
    _, org_qc_info = load_qc_table(f'{prefix.lower()}_original_failed')
    _, rep_qc_info = load_qc_table(f'{prefix.lower()}_replaced_failed')
    # Extract the unique positions (accounts for duplicated records).
    org_qc_failed_pos = np.unique(org_qc_info['pos'])
    rep_qc_failed_pos = np.unique(rep_qc_info['pos'])
    # Intialize matricies to store the effective sequence lengths.
    exons_esl = np.empty((exons_df.shape[0], 2), dtype=int)
    regions_esl = np.empty((regions_df.shape[0], 2), dtype=int)

    # Iterate through the exons.
    for i, row in exons_df.iterrows():
        # Extracts the coordinates and length.
        start, end, length = row.start, row.end, row.length
        # Intialize masks.
        org_n_failed = ((start <= org_qc_failed_pos) & (org_qc_failed_pos <= end)).sum()
        rep_n_failed = ((start <= rep_qc_failed_pos) & (rep_qc_failed_pos <= end)).sum()
        # Update the matrix.
        exons_esl[i, :] = np.array([length - org_n_failed, length - rep_n_failed])

    # Iterate through the regions.
    for i, row in regions_df.iterrows():
        # Extracts the coordinates and length.
        start, end, length = row.start, row.end, row.length
        # Intialize masks.
        org_n_failed = ((start <= org_qc_failed_pos) & (org_qc_failed_pos <= end)).sum()
        rep_n_failed = ((start <= rep_qc_failed_pos) & (rep_qc_failed_pos <= end)).sum()
        # Update the matrix.
        regions_esl[i, :] = np.array([length - org_n_failed, length - rep_n_failed])

    # Update the dataframes.
    exons_df['original_eff_seq_len'] = exons_esl[:, 0]
    exons_df['replaced_eff_seq_len'] = exons_esl[:, 1]
    regions_df['original_eff_seq_len'] = regions_esl[:, 0]
    regions_df['replaced_eff_seq_len'] = regions_esl[:, 1]
    # Export.
    exons_df.to_csv(f'./meta_data/{prefix.lower()}_ncbi_refseq_select_chr6_exons_eff_seq_lens.csv.gz', index=False)
    regions_df.to_csv(f'./meta_data/{prefix.lower()}_chr6_focal_regions_eff_seq_lens.csv.gz', index=False)
    return

# Compute the effective sequence lengths.
compute_eff_seq_lens(str(sys.argv[1]))