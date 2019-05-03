import os
from io import StringIO
import pandas as pd
from pybedtools import BedTool
from util_path import get_path


def make_snp_start_BED(snp_bed_fn):
    snp_bed = pd.read_csv(snp_bed_fn, sep="\t", header=None, usecols=[0,1,2,3], 
                          dtype={"chromStart": int, "chromEnd": int},
                          names=["chrom", "chromStart", "chromEnd", "name"])
    snp_bed.loc[:, "chromEnd"] = snp_bed.loc[:, "chromEnd"] - 1

    # BedTools require sorted BED inputs
    return snp_bed.sort_values(by=['chrom', 'chromStart'], ascending=True)

def make_TSS_start_BED(tss_bed_fn):
    tss_bed = pd.read_csv(tss_bed_fn, sep="\t", header=None, usecols=[0,1,2,3,4,5],
                          dtype={"chromStart": int, "chromEnd": int, "strand": int}, 
                          names=["chrom", "chromStart", "chromEnd", "name", "id", "strand"])

    # Determine TSS by strands
    fwd_strand = (tss_bed.loc[:, "strand"] == 1)
    rev_strand = (tss_bed.loc[:, "strand"] == -1)

    tss_bed.loc[fwd_strand, "chromEnd"] = tss_bed.loc[fwd_strand, "chromStart"]
    tss_bed.loc[rev_strand, "chromStart"] = tss_bed.loc[rev_strand, "chromEnd"]

    tss_bed = tss_bed.sort_values(by=['chrom', 'chromStart'], ascending=True)

    # BedTools require sorted BED inputs
    return tss_bed.drop_duplicates().sort_values(by=['chrom', 'chromStart'], ascending=True)

def convert_to_bedtool_obj(bed_df):
    return BedTool(bed_df.to_string(index=False, header=False, index_names=False), from_string=True)

def convert_to_df(closest_obj):
    return pd.read_csv(StringIO(str(closest_obj)), sep="\t", header=None, 
                       names=["snpChrom", "snpChromStart", "snpChromEnd", "snpName", 
                              "geneChrom", "geneChromStart", "geneChromEnd", "geneName", "geneID", "geneStrand", "distance"])


if __name__ == "__main__":
    snp_dir = get_path("vertex/SNP")
    snp_fn = "osu18_SNP.bed"
    snp_start_bed = make_snp_start_BED(os.path.join(snp_dir, snp_fn))

    gene_dir = get_path("vertex/gene")
    ensembl_TSS_fn = "Ensembl_TSS.bed"
    ensembl_TSS_start_bed = make_TSS_start_BED(os.path.join(gene_dir, ensembl_TSS_fn))
    entrez_TSS_fn = "Entrez_TSS.bed"
    entrez_TSS_start_bed = make_TSS_start_BED(os.path.join(gene_dir, entrez_TSS_fn))
    
    snp_start_bt = convert_to_bedtool_obj(snp_start_bed)
    ensembl_TSS_start_bt = convert_to_bedtool_obj(ensembl_TSS_start_bed)
    entrez_TSS_start_bt = convert_to_bedtool_obj(entrez_TSS_start_bed)

    snp_ensembl_closest = snp_start_bt.closest(ensembl_TSS_start_bt, d=True, D="b")
    snp_entrez_closest = snp_start_bt.closest(entrez_TSS_start_bt, d=True, D="b")
    
    snp_ensembl_map = convert_to_df(snp_ensembl_closest)
    snp_entrez_map = convert_to_df(snp_entrez_closest)

    res_dir = get_path("resource/NearestGene")
    snp_ensembl_map.to_csv(os.path.join(res_dir, "p1_SNP_x_NearestEnsembl.tsv"), sep="\t", header=True, index=False)
    snp_entrez_map.to_csv(os.path.join(res_dir, "p1_SNP_x_NearestEntrez.tsv"), sep="\t", header=True, index=False)

    # For any closest SNP-gene pair (s,g), if g's Ensembl ID can be mapped to an Entrez ID, then output an edge (s, g_entrez)
    # if not, output 2 edges in the edgelist file, one (s, g_ensembl) and the other (s, g_entrez)
    e2e_df = pd.read_csv(os.path.join(gene_dir, "Ensembl_x_Entrez.tsv"), sep="\t")
    
    snp_ensembl_el = snp_ensembl_map.loc[:, ["snpName", "geneID"]]
    snp_ensembl_el = snp_ensembl_el.loc[~snp_ensembl_el.geneID.isin(set(e2e_df.Ensembl_Gene_ID)), :]  # Keep Non-mappable Ensembl IDs only
    snp_entrez_el = snp_entrez_map.loc[:, ["snpName", "geneID"]]
    snp_gene_el = pd.concat([snp_ensembl_el, snp_entrez_el])

    output_dir = get_path("edge/snp-gene")
    snp_gene_el.to_csv(os.path.join(output_dir, "SNP_x_NearestGene.edgelist"), sep="\t", header=False, index=False)
    