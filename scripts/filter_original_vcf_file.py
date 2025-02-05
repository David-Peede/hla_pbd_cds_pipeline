# Import packages.
import gzip
import sys

### sys.argv[1]: path to the passed qc table with the suffix txt.gz ###
### sys.argv[2]: path to the vcf.gz file ###


# Define a function to filter the original vcf file.
def filter_original_vcf(qced_file, vcf, buffer_size=100_000):
    """
    Outputs the filtered VCF file to stdout.
    """
    # Intialize a set to store the qc'ed postions.
    qced_pos = set()
    # Open the qc'ed file.
    with gzip.open(qced_file, 'rt') as qced_data:
        # Skip the header line.
        next(qced_data)
        # For every line in the qc'ed file.
        for line in qced_data:
            # Update the qc'ed positions set.
            qced_pos.add(int(line.split()[1]))
            
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

# Filter the original vcf file.
filter_original_vcf(str(sys.argv[1]), str(sys.argv[2]))