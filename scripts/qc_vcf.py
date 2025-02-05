# Import packages.
import gzip
import sys

### sys.argv[1]: path to the qc table with the suffix txt.gz ###
### sys.argv[2]: path to the vcf.gz file ###


# Define a function to identify what duplicated entires from a single chromosome vcf.
def identify_dups(vcf):
    """
    Returns a set of duplicated positions.
    """
    # Intilailize the previous  position.
    prev_pos = None
    # Inialize a set to store duplicated positions.
    dup_set = set()
    # Read the vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line does not contains meta information.
            if not line.startswith('#'):
                # Split the line by tabs.
                spline = line.split()
                # Grab the current poistion.
                pos = spline[1]
                # If the current position is already known to be duplicated.
                if pos in dup_set:
                    # Update the previous position.
                    prev_pos = pos
                # Else-if the current position is the same as the previous position.
                elif pos == prev_pos:
                    # Add the position to the duplicated set.
                    dup_set.add(pos)
                # Else, the current position is unique.
                else:
                    # Update the previous position.
                    prev_pos = pos
    return dup_set

# Define a function to qc a line from the unfiltered vcf.
def qc_vcf_line(line, dup_set):
    """
    Returns the qc information for an autosomal TGP VCF line.
    
    QC IDs:
    - 0 = passed qc (ie biallelic snp)
    - 1 = duplicate record
    - 2 = multi alleleic
    """
    # Split the line by tabs.
    spline = line.split()
    # Grab the chromosome, position, refernce, and alternative alleles.
    chrom, pos, ref, alt = [spline[i] for i in [0, 1, 3, 4]]
    
    # If the current position has duplicated records.
    if pos in dup_set:
        # Return the qc information.
        return f'{chrom}\t{pos}\t{ref}\t{alt}\t1\n'
    # Else-if the length of refernce and alternative fields is not exactly 2,
    # ie if the site is a not mono- or biallelic snp.
    elif (len(ref) + len(alt)) != 2:
        # Return the qc information.
        return f'{chrom}\t{pos}\t{ref}\t{alt}\t2\n'
    # Else, the site is a biallelic snp.
    else:
        # Return the qc information.
        return f'{chrom}\t{pos}\t{ref}\t{alt}\t0\n'

# Define a function to parse an a single chromosome vcf file and output a qc table.
def create_qc_table(qc_table, vcf, buffer_size=100_000):
    """
    Creates a txt.gz file with the qc information per site.
    """
    # Determine the duplicated positions.
    dup_set = identify_dups(vcf)
    # Intialize a list to store the the qc information.
    qc_lines = ['chrom\tpos\tref\talt\tqc_id\n']
    
    # Read the vcf file and intilialize the qc table.
    with gzip.open(vcf, 'rt') as vcf_data, \
         gzip.open(qc_table, 'wt') as qc_file:
        # Iterate through every line in the original vcf file.
        for line in vcf_data:
            # If the line does not contains meta information.
            if not line.startswith('#'):
                # Process and append the current line.
                qc_lines.append(qc_vcf_line(line, dup_set))
            # If the number of qc lines exceeds the buffer size.
            if len(qc_lines) >= buffer_size:
                # Write the qc lines to the file.
                qc_file.writelines(qc_lines)
                # Clear the written qc lines.
                qc_lines.clear()

        # If there are still qc lines to be written.
        if qc_lines:
            # Write the remaining qc lines.
            qc_file.writelines(qc_lines)
    return

# QC the vcf file.
create_qc_table(str(sys.argv[1]), str(sys.argv[2]))