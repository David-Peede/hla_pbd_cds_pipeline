# 1000 Genomes Project (TGP) *HLA* PBD CDS Pipeline
 
 
## Notes 
 - The `vcfs` directory can be found [here](https://www.dropbox.com/scl/fo/c9mrui22qlqj08z95uhi1/AE3B3E2svflwwnjS_ATDvdc?rlkey=yqyorm09bjjst6pj3hscohxl3&st=my1k5c2i&dl=0) and to extract the `tar.xz` file run `tar -xJvf vcfs.tar.xz`.
 - I used the following packages: [BCFtoools v1.13](https://samtools.github.io/bcftools/bcftools.html), [Tabix v0.2.6](https://www.htslib.org/doc/tabix.html), [NumPy v1.22.3](https://numpy.org/doc/stable/reference/index.html), and [pandas v1.4.2](https://pandas.pydata.org/docs/).


## Code

The code was executed in the following order.


### Original TGP VCF

**‌Index the unfiltered TGP chromosome 6 VCF file.**
```bash
tabix -p vcf ./vcfs/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

**‌QC the original TGP chromosome 6 VCF file.**
```bash
python ./scripts/qc_vcf.py ./meta_data/tgp_original_qc_chr6.txt.gz ./vcfs/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```

**Extract the regions files for the original TGP chromosome 6 VCF file.‌**
```bash
python ./scripts/extract_qc_regions_files.py tgp
```

**‌Filter the original TGP chromosome 6 VCF file to only include bi-allelic SNPs.**
```bash
python ./scripts/filter_original_vcf_file.py ./meta_data/tgp_original_passed_qc_chr6.txt.gz ./vcfs/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bgzip > ./vcfs/tgp_original_chr6_no_dups_biallelic_snps.vcf.gz
```

**Index the original TGP chromosome 6 filtered VCF file.‌**
```bash
tabix -p vcf ./vcfs/tgp_original_chr6_no_dups_biallelic_snps.vcf.gz
```

**‌Subset the original TGP chromosome 6 filtered VCF file for exons only.**
```bash
bcftools view -R ./meta_data/ncbi_refseq_select_chr6_exons.txt -Oz -o ./vcfs/tgp_original_chr6_no_dups_biallelic_snps_exons.vcf.gz ./vcfs/tgp_original_chr6_no_dups_biallelic_snps.vcf.gz
```

**‌Index the original TGP chromosome 6 exons filtered VCF file.**
```bash
tabix -p vcf ./vcfs/tgp_original_chr6_no_dups_biallelic_snps_exons.vcf.gz
```


### Replaced TGP VCF

**Compile the *HLA* PBD CDS VCFs for the TGP.**
```bash
python ./scripts/compile_hla_pbd_cds_vcfs.py tgp
```

**Compile the *HLA* PBD CD QC tables for the TGP.‌**
```bash
python ./scripts/compile_hla_pbd_cds_qc_tables.py tgp
```

**‌Filter the replaced TGP chromosome 6 VCF file to only include bi-allelic SNPs.**
```bash
python ./scripts/filter_replaced_vcf_file.py ./meta_data/tgp_replaced_passed_qc_chr6.txt.gz ./meta_data/tgp_mono_bi_allelic_sites_hla_pbd_cds_vcf_no_header.txt.gz ./vcfs/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bgzip > ./vcfs/tgp_replaced_chr6_no_dups_biallelic_snps.vcf.gz
```

**‌Index the replaced TGP chromosome 6 filtered VCF file.**
```bash
tabix -p vcf ./vcfs/tgp_replaced_chr6_no_dups_biallelic_snps.vcf.gz
```

**Subset the replaced TGP chromosome 6 filtered VCF file for exons only.‌**
```bash
bcftools view -R ./meta_data/ncbi_refseq_select_chr6_exons.txt -Oz -o ./vcfs/tgp_replaced_chr6_no_dups_biallelic_snps_exons.vcf.gz ./vcfs/tgp_replaced_chr6_no_dups_biallelic_snps.vcf.gz
```

**Index the replaced TGP chromosome 6 exons filtered VCF file.**
```bash
tabix -p vcf ./vcfs/tgp_replaced_chr6_no_dups_biallelic_snps_exons.vcf.gz
```


### Subsetting VCFs

**Subset the regions for the original and replaced TGP chromosome 6 filtered VCF file.**
```bash
for dataset in original replaced; do
# HLA region.
bcftools view --regions 6:28477797-33448354 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_region.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA genes.
bcftools view --regions 6:29910309-29913647,6:31321652-31324956,6:31236526-31239869,6:32546552-32557625,6:32627244-32634434 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_genes.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# PBD exons.
bcftools view --regions 6:29910534-29910803,6:29911045-29911320,6:31323944-31324219,6:31324465-31324734,6:31238850-31239125,6:31239376-31239645,6:32551886-32552155,6:32632575-32632844 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_pbd_exons.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-A gene.
bcftools view --regions 6:29910309-29913647 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_a_gene.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-B gene.
bcftools view --regions 6:31321652-31324956 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_b_gene.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-C gene.
bcftools view --regions 6:31236526-31239869 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_c_gene.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-DRB1 gene.
bcftools view --regions 6:32546552-32557625 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_drb1_gene.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-DQB1 gene.
bcftools view --regions 6:32627244-32634434 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_dqb1_gene.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-A PBD.
bcftools view --regions 6:29910534-29910803,6:29911045-29911320 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_a_pbd_exons.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-B PBD.
bcftools view --regions 6:31323944-31324219,6:31324465-31324734 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_b_pbd_exons.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-C PBD.
bcftools view --regions 6:31238850-31239125,6:31239376-31239645 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_c_pbd_exons.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-DRB1 PBD.
bcftools view --regions 6:32551886-32552155 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_drb1_pbd_exon.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
# HLA-DQB1 PBD.
bcftools view --regions 6:32632575-32632844 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_dqb1_pbd_exon.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz
done
```

**Index all of the subsetted VCF files.‌**
```bash
for dataset in original replaced; do for subset in chr6_no_dups_biallelic_snps_hla_region chr6_no_dups_biallelic_snps_hla_genes chr6_no_dups_biallelic_snps_pbd_exons chr6_no_dups_biallelic_snps_hla_a_gene chr6_no_dups_biallelic_snps_hla_b_gene chr6_no_dups_biallelic_snps_hla_c_gene chr6_no_dups_biallelic_snps_hla_drb1_gene chr6_no_dups_biallelic_snps_hla_dqb1_gene chr6_no_dups_biallelic_snps_hla_a_pbd_exons chr6_no_dups_biallelic_snps_hla_b_pbd_exons chr6_no_dups_biallelic_snps_hla_c_pbd_exons chr6_no_dups_biallelic_snps_hla_drb1_pbd_exon chr6_no_dups_biallelic_snps_hla_dqb1_pbd_exon; do
tabix -p vcf ./vcfs/tgp_${dataset}_${subset}.vcf.gz
done; done
```

**‌Find the set difference between regions for both datasets.**
```bash
for dataset in original replaced; do
# Chromosome 6 exons \ HLA region.
bcftools isec -C ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_region.vcf.gz -n =1 -w 1 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_exons_setdiff_hla_region.vcf.gz
# HLA region \ HLA genes.
bcftools isec -C ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_region.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_genes.vcf.gz -n =1 -w 1 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_region_setdiff_hla_genes.vcf.gz
# HLA genes \ PBD exons.
bcftools isec -C ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_genes.vcf.gz ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_pbd_exons.vcf.gz -n =1 -w 1 -Oz -o ./vcfs/tgp_${dataset}_chr6_no_dups_biallelic_snps_hla_genes_setdiff_pbd_exons.vcf.gz
done
```

**Index all of the set difference VCF files.**
```bash
for dataset in original replaced; do for subset in chr6_no_dups_biallelic_snps_exons_setdiff_hla_region chr6_no_dups_biallelic_snps_hla_region_setdiff_hla_genes chr6_no_dups_biallelic_snps_hla_genes_setdiff_pbd_exons; do
tabix -p vcf ./vcfs/tgp_${dataset}_${subset}.vcf.gz
done; done
```


**Compute the effective sequence lengths per exon and focal region.**
```bash
python ./scripts/compute_effective_sequence_lengths.py tgp
```