"""
The following 2 BED files, 

- InteractorA_SNP_intxn.bed and 
- InteractorA_SNP_intxn.bed 

have 9 columns each (without a header row):

| 0              | 1                | 2              | 3              | 4       | 5         | 6         | 7          | 8       |
| Interactor_Chr | Interactor_Start | Interactor_End | Interaction_ID | SNP_Chr | SNP_Start | SNP_End   | SNP_ID     | Overlap |
|----------------|------------------|----------------|----------------|---------|-----------|-----------|------------|---------|
| chr7           | 155604229        | 155606229      | 18695675_1     | chr7    | 155604940 | 155604941 | rs9333594  | 1       |
| chr20          | 57456862         | 57471567       | 24413736_155   | chr20   | 57458103  | 57458104  | rs60865276 | 1       |
| chr20          | 57456862         | 57471567       | 24413736_155   | chr20   | 57465570  | 57465571  | rs6123837  | 1       |

==========================

The following 2 BED files:

- InteractorA_Ensembl_TSS_intxn.bed,
- InteractorA_Ensembl_TSS_intxn.bed

have 11 columns (without a header row):

| 0              | 1                | 2              | 3              | 4        | 5         | 6       | 7         | 8               | 9      | 10      |
| Interactor_Chr | Interactor_Start | Interactor_End | Interaction_ID | Gene_Chr | TSS_Start | TSS_End | Gene_Name | Gene_Ensembl_ID | Strand | Overlap |
|----------------|------------------|----------------|----------------|----------|-----------|---------|-----------|-----------------|--------|---------|
| chr11          | 1783071          | 1788116        | 19890323_1     | chr11    | 1785221   | 1785222 | CTSD      | ENSG00000117984 | -1     | 1       |
| chr11          | 1783071          | 1788116        | 19890323_1     | chr11    | 1785130   | 1785131 | CTSD      | ENSG00000117984 | -1     | 1       |
| chr11          | 1783071          | 1788116        | 19890323_1     | chr11    | 1783632   | 1783633 | CTSD      | ENSG00000117984 | -1     | 1       |

==========================

The following 2 BED files:

- InteractorA_Ensembl_promoter_intxn.bed,
- InteractorA_Ensembl_promoter_intxn.bed

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

def get_genes_with_interacted_TSSs(interactor_tss_intxn_bed):
    """
    Read the `BedTools intersect` output of interactor and TSS BEDs, output IDs of genes whose TSS overlaps with an interator region
    """
    return pd.read_csv(interactor_tss_intxn_bed, header=None, sep="\t", usecols=[3, 7, 8],
                       names=["interaction_id", "gene_name", "gene_ensembl_id"]).drop_duplicates()

def get_genes_with_interacted_promoters(intx_prm_intersect_bed):
    """
    Read the `BedTools intersect` output of interactor and promoter BEDs, output IDs of genes whose promoter overlaps with an interator region
    """
    return get_genes_with_interacted_TSSs(intx_prm_intersect_bed)  # delegation


def get_interacted_SNPs(interactor_snp_intxn_bed):
    """
    Read the `BedTools intersect` output of interactor and SNP BEDs, output IDs of SNPs who overlaps with an interator region
    """
    return pd.read_csv(interactor_snp_intxn_bed, header=None, sep="\t", usecols=[3, 7],
                       names=["interaction_id", "rs_id"]).drop_duplicates()

def _cross_inner_join(xA_df, xB_df, yA_df, yB_df):
    """
    x, y: the types of entities that reside in 4DGenome interactor regions (i.e. "SNP", "TSS" and "Promoter")
    A, B: interactor regions

    E.g. if x = "SNP", then "xA" means SNPs that reside in interactor A regions

    Perform 2 joins on "interaction_id", one between (xA, yB) and the other between (xB, yA); 
    then concatenate the results of these two joins

    xA ==+ +>> yA
          X
    xB >>+ +== yB
    """
    xA_yB_map = xA_df.merge(yB_df, on="interaction_id", how="inner")
    xB_yA_map = xB_df.merge(yA_df, on="interaction_id", how="inner")

    xy_map = pd.concat([xA_yB_map, xB_yA_map], axis=0, ignore_index=True).drop_duplicates()

    return xy_map

def _inner_join(xA_df, xB_df):
    """
    Simple innner join on "interaction_id"
    """
    xA_xB_map = xA_df.merge(xB_df, on="interaction_id", how="inner", suffixes=('_A', '_B'))

    return xA_xB_map

def output_snp_gene_edgelist(snp_gene_map):             
    return snp_gene_map.loc[:, ["rs_id", "gene_ensembl_id"]].drop_duplicates()

def output_snp_snp_edgelist(snp_snp_map):
    return snp_snp_map.loc[:, ["rs_id_A", "rs_id_B"]].drop_duplicates()


if __name__ == "__main__":
    _4DGenome_dir = get_path("resource/4DGenome")

    snp_in_interactorA = get_interacted_SNPs(os.path.join(_4DGenome_dir, "InteractorA_SNP_intxn.bed"))
    snp_in_interactorB = get_interacted_SNPs(os.path.join(_4DGenome_dir, "InteractorB_SNP_intxn.bed"))

    gene_tss_in_interactorA = get_genes_with_interacted_TSSs(os.path.join(_4DGenome_dir, "InteractorA_Ensembl_TSS_intxn.bed"))
    gene_tss_in_interactorB = get_genes_with_interacted_TSSs(os.path.join(_4DGenome_dir, "InteractorB_Ensembl_TSS_intxn.bed"))

    gene_prm_in_interactorA = get_genes_with_interacted_promoters(os.path.join(_4DGenome_dir, "InteractorA_Ensembl_promoter_intxn.bed"))
    gene_prm_in_interactorB = get_genes_with_interacted_promoters(os.path.join(_4DGenome_dir, "InteractorB_Ensembl_promoter_intxn.bed"))

    # interacted (SNP, gene) pairs (via TSSs/promoters)
    snp_tss_map = _cross_inner_join(snp_in_interactorA, snp_in_interactorB, gene_tss_in_interactorA, gene_tss_in_interactorB)
    snp_prm_map = _cross_inner_join(snp_in_interactorA, snp_in_interactorB, gene_prm_in_interactorA, gene_prm_in_interactorB)
    # interacted (SNP, SNP) pairs
    snp_snp_map = _inner_join(snp_in_interactorA, snp_in_interactorB)
    
    snp_tss_map.to_csv(os.path.join(_4DGenome_dir, "p1_SNP_x_4DGenome_TSS.tsv"), sep="\t", index=False)
    snp_prm_map.to_csv(os.path.join(_4DGenome_dir, "p1_SNP_x_4DGenome_promoter.tsv"), sep="\t", index=False)
    snp_snp_map.to_csv(os.path.join(_4DGenome_dir, "p1_SNP_x_4DGenome_SNP.tsv"), sep="\t", index=False)

    snp_tss_el = output_snp_gene_edgelist(snp_tss_map)
    snp_prm_el = output_snp_gene_edgelist(snp_prm_map)
    snp_snp_el = output_snp_snp_edgelist(snp_snp_map)
    
    _output_dir = get_path("edge/snp-gene")
    snp_tss_el.to_csv(os.path.join(_output_dir, "SNP_x_4DGenome_TSS.edgelist"), sep="\t", index=False, header=False)
    snp_prm_el.to_csv(os.path.join(_output_dir, "SNP_x_4DGenome_promoter.edgelist"), sep="\t", index=False, header=False)
    
    _output_dir = get_path("edge/snp-snp")
    snp_snp_el.to_csv(os.path.join(_output_dir, "SNP_x_4DGenome_SNP.edgelist"), sep="\t", index=False, header=False)