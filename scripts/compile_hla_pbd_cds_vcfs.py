# Import packages.
import numpy as np
import pandas as pd
import sys

### sys.argv[1]: dataset prefix ###


# Define a function to load the targt tables.
def load_targt_table(prefix, suffix):
    """
    Returns the HLA PBD CDS table as a pandas dataframe.
    """
    return pd.read_csv(f'./hla_pbd_cds_tables/{prefix.upper()}full_HLASeq{suffix.upper()}.csv')

# Define a function to compile and qc the targt tables.
def compile_hla_pbd_cds_tables(prefix):
    """
    Returns headerless vcf files as a .txt.gz .
    """
    # Load and sort all of the hla tables.
    hla_df = pd.concat(
        [load_targt_table(prefix, locus) for locus in ['a', 'b', 'c', 'dqb1', 'drb1']]
    ).sort_values(by=['POS'], ascending=True, ignore_index=True)
    
    # Subset the dataframe to export later.
    vcf_df = hla_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']].copy()
    # Reformat the info column.
    vcf_df['INFO'] = 'VT=' + vcf_df['INFO']
    
    # Extract the genotypes as a numpy matrix.
    genos = hla_df.to_numpy()[:, 9:]
    # Extract the sample ids.
    samp_ids = hla_df.columns.values[9:]
    # Extract the reference alleles.
    refs = hla_df.REF.values
    # Intialize a list to store the alternative alleles.
    alts = []
    
    # For every position.
    for i, ref in enumerate(refs):
        # Determine the unique alleles.
        alleles = np.unique(genos[i, :])
        # Recode the reference alleles.
        genos[i, :][genos[i, :] == ref] = '0'
        # If there are missing genotypes.
        if 'X' in alleles:
            # Recode the missing alleles.
            genos[i, :][genos[i, :] == 'X'] = '.'
        # Remove the refernce and missing alleles.
        alleles = alleles[(alleles != ref) & (alleles != 'X')]
        # If there are no alleles.
        if alleles.size == 0:
            # Update the alternate allele column.
            alts.append('.')
        # Else.
        else:
            # For every alternative allele.
            for j, alt in enumerate(alleles):
                # Recode the alternative allele.
                genos[i, :][genos[i, :] == alt] = f'{int(j+1)}'
            # Update the alternate allele column.
            alts.append(','.join(alleles))
    
    # Update the alternative allele column.
    vcf_df['ALT'] = alts
    # Create a mask for the mono- and bi-allelic sites.
    is_qced = np.array([len(alt) == 1 for alt in alts])
    
    # For every sample id.
    for i, samp_id in enumerate(samp_ids):
        # If this is the first chromsome.
        if '.' not in samp_id:
            # Update the vcf dataframe.
            vcf_df[samp_id] = [f'{hap1}/{hap2}' for hap1, hap2 in zip(genos[:, i], genos[:, i+1])]
    
    # Export the vcf dataframes.
    vcf_df.to_csv(f'./meta_data/{prefix.lower()}_all_sites_hla_pbd_cds_vcf_no_header.txt.gz', sep='\t', index=False)
    vcf_df[is_qced].to_csv(f'./meta_data/{prefix.lower()}_mono_bi_allelic_sites_hla_pbd_cds_vcf_no_header.txt.gz', sep='\t', index=False)
    vcf_df[~is_qced].to_csv(f'./meta_data/{prefix.lower()}_multi_allelic_sites_hla_pbd_cds_vcf_no_header.txt.gz', sep='\t', index=False)
    return

# Compile the TARGT tables.
compile_hla_pbd_cds_tables(str(sys.argv[1]))