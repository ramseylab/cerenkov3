## GTEx Download

- [GTEx Analysis V7 (dbGaP Accession phs000424.v7.p2)](https://gtexportal.org/home/datasets#datasetDiv1)
    - => [Single-Tissue cis-eQTL Data](https://gtexportal.org/home/datasets#filesetFilesDiv15)
        - => eGene and significant variant-gene associations based on permutations: [GTEx_Analysis_v7_eQTL.tar.gz](https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz)

## GTEx File Headers

See https://github.com/intermine/intermine/issues/1607

### All SNP-gene pairs (File extension: `*.egenes.txt.gz`)

Note that the `*.egenes.txt.gz` files contain data for all genes tested; to obtain the list of eGenes, select the rows with `qval <= 0.05`.
  
Column headers in the file:

| Column Name              | Comment                                                                                       | 
|--------------------------|-----------------------------------------------------------------------------------------------| 
| gene_id                  | GENCODE/Ensembl gene ID                                                                       | 
| gene_name                | GENCODE gene name                                                                             | 
| gene_chr                 | chromosome (gene)                                                                             | 
| gene_start               | gene start position (in base pairs; 1-based coordinates)                                      | 
| gene_end                 | gene end position (in base pairs; 1-based coordinates)                                        | 
| strand                   | genomic strand                                                                                | 
| num_var                  | number of variants in cis-window                                                              | 
| beta_shape1              | 1st shape parameter of the fitted Beta distribution: B(shape1, shape2)                        | 
| beta_shape2              | 2nd shape parameter of the fitted Beta distribution: B(shape1, shape2)                        | 
| true_df                  | Effective degrees of freedom the Beta distribution approximation                              | 
| variant_id               | variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37                   | 
| tss_distance             | Distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise | 
| chr                      | chromosome (variant; same as gene_chr for cis-eQTLs)                                          | 
| snp_pos                  | position of the first reference base of the variant                                           | 
| ref                      | reference sequence of the variant                                                             | 
| alt                      | alternate sequence of the variant                                                             | 
| rs_id_dbSNP142_GRCh37p13 | dbSNP142 rsID                                                                                 | 
| num_alt_per_site         | number of alternative alleles observed at this site                                           | 
| minor_allele_samples     | number of samples carrying the minor allele                                                   | 
| minor_allele_count       | total number of minor alleles across individuals                                              | 
| maf                      | minor allele frequency observed in the set of donors for a given tissue                       | 
| ref_factor               | '1', when the minor allele is the alt base, '-1' when the minor allele is the reference base  | 
| pval_nominal             | nominal p-value associated with the most significant variant for this gene                    | 
| slope                    | regression slope                                                                              | 
| slope_se                 | standard error of the regression slope                                                        | 
| pval_perm                | permutation p-value                                                                           | 
| pval_beta                | beta-approximated permutation p-value                                                         | 
| qval                     | Storey q-value derived from pval_beta                                                         | 
| pval_nominal_threshold   | nominal p-value threshold for calling a variant-gene pair significant for the gene            | 

### Significant SNP-gene pairs (File extension: *.signif_snpgene_pairs.txt.gz)

Column headers in the file:

| Column Name            | Comment                                                                                       | 
|------------------------|-----------------------------------------------------------------------------------------------| 
| variant_id             | variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37                   | 
| gene_id                | GENCODE/Ensembl gene ID                                                                       | 
| tss_distance           | distance between variant and transcription start site. Positive when variant is downstream of the TSS, negative otherwise | 
| pval_nominal           | nominal p-value                                                                               | 
| slope                  | regression slope                                                                              | 
| slope_se               | standard error of the regression slope                                                        | 
| slope_fpkm             | regression slope in FPKM units, calculated using quantile normalized expression tables        | 
| slope_fpkm_se          | standard error of slope_fpkm                                                                  | 
| pval_nominal_threshold | nominal p-value threshold for calling a variant-gene pair significant for the gene            | 
| min_pval_nominal       | smallest nominal p-value for the gene                                                         | 
| pval_beta              | beta-approximated permutation p-value for the gene                                            | 

## Data Processing

Here we will simply match egenes to our known SNPs and output SNP-gene edges.