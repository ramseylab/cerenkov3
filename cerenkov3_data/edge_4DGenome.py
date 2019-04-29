"""
The following 2 BED files, 

- 4DGenome_InteractorA_SNP.bed and 
- 4DGenome_InteractorB_SNP.bed 

have 9 columns each (without a header row):

| 0              | 1                | 2              | 3              | 4       | 5         | 6         | 7          | 8       |
| Interactor_Chr | Interactor_Start | Interactor_End | Interaction_ID | SNP_Chr | SNP_Start | SNP_End   | SNP_ID     | Overlap |
|----------------|------------------|----------------|----------------|---------|-----------|-----------|------------|---------|
| chr7           | 155604229        | 155606229      | 18695675_1     | chr7    | 155604940 | 155604941 | rs9333594  | 1       |
| chr20          | 57456862         | 57471567       | 24413736_155   | chr20   | 57458103  | 57458104  | rs60865276 | 1       |
| chr20          | 57456862         | 57471567       | 24413736_155   | chr20   | 57465570  | 57465571  | rs6123837  | 1       |

==========================

The following 2 BED files:

- 4DGenome_InteractorA_Ensembl_TSS.bed,
- 4DGenome_InteractorB_Ensembl_TSS.bed

have 11 columns (without a header row):

| 0              | 1                | 2              | 3              | 4        | 5         | 6       | 7         | 8               | 9      | 10      |
| Interactor_Chr | Interactor_Start | Interactor_End | Interaction_ID | Gene_Chr | TSS_Start | TSS_End | Gene_Name | Gene_Ensembl_ID | Strand | Overlap |
|----------------|------------------|----------------|----------------|----------|-----------|---------|-----------|-----------------|--------|---------|
| chr11          | 1783071          | 1788116        | 19890323_1     | chr11    | 1785221   | 1785222 | CTSD      | ENSG00000117984 | -1     | 1       |
| chr11          | 1783071          | 1788116        | 19890323_1     | chr11    | 1785130   | 1785131 | CTSD      | ENSG00000117984 | -1     | 1       |
| chr11          | 1783071          | 1788116        | 19890323_1     | chr11    | 1783632   | 1783633 | CTSD      | ENSG00000117984 | -1     | 1       |

==========================

The following 2 BED files:

- 4DGenome_InteractorA_Ensembl_promoter.bed,
- 4DGenome_InteractorB_Ensembl_promoter.bed

have 11 columns (without a header row):

| 0              | 1                | 2              | 3              | 4        | 5              | 6            | 7         | 8               | 10     | 11      |
| Interactor_Chr | Interactor_Start | Interactor_End | Interaction_ID | Gene_Chr | Promoter_Start | Promoter_End | Gene_Name | Gene_Ensembl_ID | Strand | Overlap |
|----------------|------------------|----------------|----------------|----------|----------------|--------------|-----------|-----------------|--------|---------|
| chr11          | 61094562         | 61096562       | 22791839_3     | chr11    | 61093200       | 61095701     | DDB1      | ENSG00000167986 | -1     | 1139    |
| chr11          | 61094562         | 61096562       | 22791839_3     | chr11    | 61094522       | 61097023     | DDB1      | ENSG00000167986 | -1     | 2000    |
| chr11          | 61094562         | 61096562       | 22791839_3     | chr11    | 61093895       | 61096396     | DDB1      | ENSG00000167986 | -1     | 1834    |
"""

import os
import pandas as pd
from util_path import get_path

def read_intx_tss_bed(filename):
    return pd.read_csv(filename, header=None, sep="\t", usecols=[3, 7, 8],
                       names=["interaction_id", "gene_name", "gene_ensembl_id"]).drop_duplicates()

def read_intx_promoter_bed(filename):
    return read_intx_tss_bed(filename)  # delegation


def read_intx_snp_bed(filename):
    return pd.read_csv(filename, header=None, sep="\t", usecols=[3, 7],
                       names=["interaction_id", "rs_id"]).drop_duplicates()

def intx_cross_inner_join(xA_df, xB_df, yA_df, yB_df):
    """
    x, y: the types of entities that reside in 4DGenome interactor regions (i.e. "SNP", "TSS" and "Promoter")
    A, B: interactor regions

    E.g. if x = "SNP", then "xA" means SNPs that reside in interactor A regions
    """
    xA_yB_map = xA_df.merge(yB_df, on="interaction_id", how="inner")
    xB_yA_map = xB_df.merge(yA_df, on="interaction_id", how="inner")

    xy_map = pd.concat([xA_yB_map, xB_yA_map], axis=0, ignore_index=True).drop_duplicates()

    return xy_map

def intx_inner_join(xA_df, xB_df):
    xA_xB_map = xA_df.merge(xB_df, on="interaction_id", how="inner", suffixes=('_A', '_B'))

    return xA_xB_map

def output_snp_gene_edgelist(snp_gene_map):             
    return snp_gene_map.loc[:, ["rs_id", "gene_ensembl_id"]].drop_duplicates()

def output_snp_snp_edgelist(snp_snp_map):
    return snp_snp_map.loc[:, ["rs_id_A", "rs_id_B"]].drop_duplicates()


if __name__ == "__main__":
    _4DGenome_dir = get_path("edge/snp-gene/Chromatin-Interaction/4DGenome")

    snp_A_df = read_intx_snp_bed(os.path.join(_4DGenome_dir, "4DGenome_InteractorA_SNP.bed"))
    snp_B_df = read_intx_snp_bed(os.path.join(_4DGenome_dir, "4DGenome_InteractorB_SNP.bed"))

    tss_A_df = read_intx_tss_bed(os.path.join(_4DGenome_dir, "4DGenome_InteractorA_Ensembl_TSS.bed"))
    tss_B_df = read_intx_tss_bed(os.path.join(_4DGenome_dir, "4DGenome_InteractorB_Ensembl_TSS.bed"))

    prm_A_df = read_intx_promoter_bed(os.path.join(_4DGenome_dir, "4DGenome_InteractorA_Ensembl_promoter.bed"))
    prm_B_df = read_intx_promoter_bed(os.path.join(_4DGenome_dir, "4DGenome_InteractorB_Ensembl_promoter.bed"))

    # SNP-gene edges (via TSSs/promoters)
    snp_tss_map = intx_cross_inner_join(snp_A_df, snp_B_df, tss_A_df, tss_B_df)
    snp_prm_map = intx_cross_inner_join(snp_A_df, snp_B_df, prm_A_df, prm_B_df)
    
    snp_tss_map.to_csv(os.path.join(_4DGenome_dir, "p1_SNP_x_4DGenome_TSS.tsv"), sep="\t", index=False)
    snp_prm_map.to_csv(os.path.join(_4DGenome_dir, "p1_SNP_x_4DGenome_promoter.tsv"), sep="\t", index=False)

    snp_tss_el = output_snp_gene_edgelist(snp_tss_map)
    snp_prm_el = output_snp_gene_edgelist(snp_prm_map)
    
    snp_tss_el.to_csv(os.path.join(_4DGenome_dir, "p2_SNP_x_4DGenome_TSS.edgelist"), sep="\t", index=False, header=False)
    snp_prm_el.to_csv(os.path.join(_4DGenome_dir, "p2_SNP_x_4DGenome_promoter.edgelist"), sep="\t", index=False, header=False)
    
    # SNP-SNP edges
    _4DGenome_dir = get_path("edge/snp-snp/Chromatin-Interaction/4DGenome")

    snp_snp_map = intx_inner_join(snp_A_df, snp_B_df)
    snp_snp_map.to_csv(os.path.join(_4DGenome_dir, "p1_SNP_x_4DGenome_SNP.tsv"), sep="\t", index=False)
    
    snp_snp_el = output_snp_snp_edgelist(snp_snp_map)
    snp_snp_el.to_csv(os.path.join(_4DGenome_dir, "p2_SNP_x_4DGenome_SNP.edgelist"), sep="\t", index=False, header=False)