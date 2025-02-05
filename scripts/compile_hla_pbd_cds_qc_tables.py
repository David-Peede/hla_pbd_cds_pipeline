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

# Define a function to load the hla pbd cds vcf files.
def load_hla_pbd_cds_vcf(prefix):
    """
    Returns the HLA PBD CDS vcf as a pandas dataframe.
    """
    return pd.read_csv(f'./meta_data/{prefix}_hla_pbd_cds_vcf_no_header.txt.gz', sep='\t')

# Define a function to generate the qc tables for the targt vcf.
def compile_hla_pbd_cds_qc_tables(prefix):
    """
    Returns the HLA PBD CDS replaced qc tables.
    """
    # Load the all sites hla pbd cds table.
    hla_pbd_cds_df = load_hla_pbd_cds_vcf(f'{prefix.lower()}_all_sites')
    # Extract the poistions, reference, and alternative allele columns.
    hla_pbd_cds_pos = hla_pbd_cds_df.POS.values
    hla_pbd_cds_refs = hla_pbd_cds_df.REF.values
    hla_pbd_cds_alts = hla_pbd_cds_df.ALT.values
    
    # Intialize allelic masks.
    is_mono = hla_pbd_cds_alts == '.'
    is_bi = np.array([len(alt) == 1 for alt in hla_pbd_cds_alts]) & ~is_mono
    is_tri = np.array([alt.count(',') == 1 for alt in hla_pbd_cds_alts])
    is_quad = np.array([alt.count(',') == 2 for alt in hla_pbd_cds_alts])
    
    # Intialize a qc id column.
    hla_pbd_cds_qc_ids = np.full(hla_pbd_cds_pos.size, -1)
    hla_pbd_cds_qc_ids[is_bi] = 0
    hla_pbd_cds_qc_ids[is_tri | is_quad] = 2
    # Intialize a dictionary.
    hla_pbd_cds_qc_dicc = {
        'chrom': np.full(hla_pbd_cds_pos.size, 6),
        'pos': hla_pbd_cds_pos,
        'ref': hla_pbd_cds_refs,
        'alt': hla_pbd_cds_alts,
        'qc_id': hla_pbd_cds_qc_ids,
    }
    
    # Load the original qc report.
    org_qc_df, org_qc_info = load_qc_table(f'{prefix.lower()}_original')
    # Determine the sites that will be replaced.
    org_is_same = np.isin(org_qc_info['pos'], hla_pbd_cds_pos)
    
    # Concatenate the qc dataframes.
    hla_pbd_cds_qc_df = pd.concat(
        [org_qc_df[~org_is_same], pd.DataFrame(hla_pbd_cds_qc_dicc)]
    ).sort_values(by=['pos'], ascending=True, ignore_index=True)
    # Intialize the qc'ed mask.
    is_qced = hla_pbd_cds_qc_df.qc_id.values <= 0
    # Export the qc results.
    hla_pbd_cds_qc_df.to_csv(f'./meta_data/{prefix.lower()}_replaced_qc_chr6.txt.gz', sep='\t', index=False)
    hla_pbd_cds_qc_df[is_qced].to_csv(f'./meta_data/{prefix.lower()}_replaced_passed_qc_chr6.txt.gz', sep='\t', index=False)
    hla_pbd_cds_qc_df[~is_qced].to_csv(f'./meta_data/{prefix.lower()}_replaced_failed_qc_chr6.txt.gz', sep='\t', index=False)
    return

# Generate the targt qc tables.
compile_hla_pbd_cds_qc_tables(str(sys.argv[1]))