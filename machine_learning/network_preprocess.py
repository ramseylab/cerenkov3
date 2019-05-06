import os
import pandas as pd


def read_snp_gene_edge_list(fn):
    print("Loading {}".format(fn))
    return pd.read_csv(fn, sep="\t", header=None, names=["rsID", "GeneID"],
                       dtype={"rsID": str, "GeneID": str}).drop_duplicates()

def read_gene_gene_edge_list(fn):
    print("Loading {}".format(fn))
    return pd.read_csv(fn, sep="\t", header=None, names=["GeneID_A", "GeneID_B"],
                       dtype={"GeneID_A": str, "GeneID_B": str}).drop_duplicates()

def read_all_snp_gene_edge_lists():
    snp_gene_edge_dir = "../cerenkov3_data/edge/snp-gene"
    
    fn_4DG_tss = "SNP_x_4DGenome_TSS.edgelist"
    fn_4DG_prm = "SNP_x_4DGenome_promoter.edgelist"
    fn_GTEx    = "SNP_x_GTEx.edgelist"
    fn_TFBS    = "SNP_x_EncodeTFBS.edgelist"
    fn_NG      = "SNP_x_NearestGene.edgelist"

    el_4DG_tss = read_snp_gene_edge_list(os.path.join(snp_gene_edge_dir, fn_4DG_tss))
    el_4DG_prm = read_snp_gene_edge_list(os.path.join(snp_gene_edge_dir, fn_4DG_prm))
    el_GTEx    = read_snp_gene_edge_list(os.path.join(snp_gene_edge_dir, fn_GTEx))
    el_TFBS    = read_snp_gene_edge_list(os.path.join(snp_gene_edge_dir, fn_TFBS))
    el_NG      = read_snp_gene_edge_list(os.path.join(snp_gene_edge_dir, fn_NG))

    return el_4DG_tss, el_4DG_prm, el_GTEx, el_TFBS, el_NG

def read_all_gene_gene_edge_lists():
    gene_gene_edge_dir = "../cerenkov3_data/edge/gene-gene"

    fn_coexpedia = "Coexpedia.edgelist"
    fn_humannet  = "HumanNet.edgelist"
    fn_biogrid   = "BioGRID.edgelist"

    el_coexpedia = read_gene_gene_edge_list(os.path.join(gene_gene_edge_dir, fn_coexpedia))
    el_humannet  = read_gene_gene_edge_list(os.path.join(gene_gene_edge_dir, fn_humannet))
    el_biogrid = read_gene_gene_edge_list(os.path.join(gene_gene_edge_dir, fn_biogrid))

    return el_coexpedia, el_humannet, el_biogrid

def get_unique_snp_gene_IDs(*sg_els):
    """
    sg_els == SNP-gene edge lists
    """
    el_composed = pd.concat(sg_els, axis=0).drop_duplicates()

    snp_id_df = el_composed.loc[:, ["rsID"]].drop_duplicates().reset_index(drop=True)

    is_ensembl = el_composed.GeneID.str.startswith("E")
    ensembl_id_df = el_composed.loc[is_ensembl, ["GeneID"]].drop_duplicates().reset_index(drop=True)
    entrez_id_df = el_composed.loc[~is_ensembl, ["GeneID"]].drop_duplicates().reset_index(drop=True)

    snp_id_df.rename(columns={"rsID": "ID"}, inplace=True)
    ensembl_id_df.rename(columns={"GeneID": "ID"}, inplace=True)
    entrez_id_df.rename(columns={"GeneID": "ID"}, inplace=True)

    return snp_id_df, ensembl_id_df, entrez_id_df

def subset_gene_gene_edge_list(gg_el, entrez_id_set):
    subset_flag = gg_el.GeneID_A.isin(entrez_id_set) & gg_el.GeneID_B.isin(entrez_id_set)

    return gg_el.loc[subset_flag, :].sort_values(by=["GeneID_A", "GeneID_B"])

def assign_int_id(*id_dfs):
    starts_from = 1
    
    for id_df in id_dfs:
        yield id_df.assign(INT_ID=id_df.index + starts_from)

        starts_from += len(id_df.index)

def apply_int_id_to_snp_gene_edge_list(sg_el, snp_id_df, gene_id_df):
    sg_el = sg_el.merge(snp_id_df, how="left", left_on="rsID", right_on="ID")
    sg_el = sg_el.merge(gene_id_df, how="left", left_on="GeneID", right_on="ID", suffixes=('_a', '_b'))

    return sg_el.loc[:, ["INT_ID_a", "INT_ID_b"]]

def apply_int_id_to_gene_gene_edge_list(gg_el, gene_id_df):
    gg_el = gg_el.merge(gene_id_df, how="left", left_on="GeneID_A", right_on="ID")
    gg_el = gg_el.merge(gene_id_df, how="left", left_on="GeneID_B", right_on="ID", suffixes=('_a', '_b'))

    return gg_el.loc[:, ["INT_ID_a", "INT_ID_b"]]

if __name__ == "__main__":
    out_dir = "./INT_ID_EDGELIST"

    el_4DG_tss, el_4DG_prm, el_GTEx, el_TFBS, el_NG = read_all_snp_gene_edge_lists()
    el_coexpedia, el_humannet, el_biogrid = read_all_gene_gene_edge_lists()

    snp_id_df, ensembl_id_df, entrez_id_df = get_unique_snp_gene_IDs(el_4DG_tss, el_4DG_prm, el_GTEx, el_TFBS, el_NG)
    snp_id_df, ensembl_id_df, entrez_id_df = assign_int_id(snp_id_df, ensembl_id_df, entrez_id_df)
    gene_id_df = pd.concat([ensembl_id_df, entrez_id_df], axis=0)

    snp_id_df.loc[:, ["ID", "INT_ID"]].to_csv(os.path.join(out_dir, "SNP_INT_ID.tsv"), sep="\t", index=False, header=True)
    ensembl_id_df.loc[:, ["ID", "INT_ID"]].to_csv(os.path.join(out_dir, "Ensembl_Gene_INT_ID.tsv"), sep="\t", index=False, header=True)
    entrez_id_df.loc[:, ["ID", "INT_ID"]].to_csv(os.path.join(out_dir, "Entrez_Gene_INT_ID.tsv"), sep="\t", index=False, header=True)

    el_4DG_tss = apply_int_id_to_snp_gene_edge_list(el_4DG_tss, snp_id_df, gene_id_df)
    el_4DG_prm = apply_int_id_to_snp_gene_edge_list(el_4DG_prm, snp_id_df, gene_id_df)
    el_GTEx    = apply_int_id_to_snp_gene_edge_list(el_GTEx, snp_id_df, gene_id_df)
    el_TFBS    = apply_int_id_to_snp_gene_edge_list(el_TFBS, snp_id_df, gene_id_df)
    el_NG      = apply_int_id_to_snp_gene_edge_list(el_NG, snp_id_df, gene_id_df)

    el_4DG_tss.to_csv(os.path.join(out_dir, "SNP_x_4DGenome_TSS.INT_ID.edgelist") , sep="\t", header=False, index=False)
    el_4DG_prm.to_csv(os.path.join(out_dir, "SNP_x_4DGenome_promoter.INT_ID.edgelist") , sep="\t", header=False, index=False)
    el_GTEx.to_csv(os.path.join(out_dir, "SNP_x_GTEx.INT_ID.edgelist") , sep="\t", header=False, index=False)
    el_TFBS.to_csv(os.path.join(out_dir, "SNP_x_EncodeTFBS.INT_ID.edgelist") , sep="\t", header=False, index=False)
    el_NG.to_csv(os.path.join(out_dir, "SNP_x_NearestGene.INT_ID.edgelist") , sep="\t", header=False, index=False)
    
    entrez_id_set = set(entrez_id_df.ID)
    el_coexpedia = subset_gene_gene_edge_list(el_coexpedia, entrez_id_set)
    el_humannet  = subset_gene_gene_edge_list(el_humannet, entrez_id_set)
    el_biogrid   = subset_gene_gene_edge_list(el_biogrid, entrez_id_set)

    el_coexpedia = apply_int_id_to_gene_gene_edge_list(el_coexpedia, gene_id_df)
    el_humannet  = apply_int_id_to_gene_gene_edge_list(el_humannet, gene_id_df)
    el_biogrid   = apply_int_id_to_gene_gene_edge_list(el_biogrid, gene_id_df)

    el_coexpedia.to_csv(os.path.join(out_dir, "Coexpedia.INT_ID.edgelist") , sep="\t", header=False, index=False)
    el_humannet.to_csv(os.path.join(out_dir, "HumanNet.INT_ID.edgelist") , sep="\t", header=False, index=False)
    el_biogrid.to_csv(os.path.join(out_dir, "BioGRID.INT_ID.edgelist") , sep="\t", header=False, index=False)
    