import os
from functools import reduce
import pandas as pd
from util_path import get_path
from util_ensembl import map_Ensembl_IDs_to_Entrez


def read_snp_df(fn):
    return pd.read_csv(fn, header=None, sep="\t", usecols=[3], names=["rs_id"])

def list_files(directory, fn_suffix):
    """
    list all files in a certain directory ending with a suffix
    """
    return [f for f in os.listdir(directory) if f.endswith(fn_suffix)]

def read_eqtl_df(fn):
    eqtl_df = pd.read_csv(fn, sep="\t", usecols=["gene_id", "gene_name", "rs_id_dbSNP147_GRCh37p13", "pval_perm", "pval_beta"])
    eqtl_df = eqtl_df.rename({"gene_id": "gene_ensembl_id",
                              "rs_id_dbSNP147_GRCh37p13": "rs_id"}, axis=1)
    return eqtl_df

def gen_snp_egene_map(snp_df, eqtl_dir, eqtl_suffix):
    eqtl_fn_list = list_files(eqtl_dir, eqtl_suffix)

    for eqtl_fn in eqtl_fn_list:
        eqtl_path = os.path.join(eqtl_dir, eqtl_fn)        
        eqtl_df = read_eqtl_df(eqtl_path)
        
        snp_egene_map = snp_df.merge(eqtl_df, on="rs_id", how="inner")
        
        eqtl_source = eqtl_fn.split(eqtl_suffix)[0]
        snp_egene_map = snp_egene_map.assign(eqtl_source = eqtl_source)
        
        yield snp_egene_map


if __name__ == "__main__":
    snp_dir = get_path("vertex/SNP")
    snp_fn = "osu18_SNP.bed"
    snp_df = read_snp_df(os.path.join(snp_dir, snp_fn))

    res_dir = get_path("resource/GTEx")
    eqtl_dir = os.path.join(res_dir, "GTEx_Analysis_v7_eQTL")
    eqtl_suffix = ".v7.egenes.txt"

    snp_egene_map = pd.concat(gen_snp_egene_map(snp_df, eqtl_dir, eqtl_suffix), ignore_index=True)
    # Get rid of version numbers
    snp_egene_map.loc[:, "gene_ensembl_id"] = snp_egene_map.loc[:, "gene_ensembl_id"].apply(lambda x: x.split(".")[0])
    
    snp_egene_map.to_csv(os.path.join(res_dir, "p1_SNP_x_GTEx.tsv"), sep="\t", index=False)

    snp_egene_el = snp_egene_map.loc[:, ["rs_id", "gene_ensembl_id"]].drop_duplicates()
    snp_egene_el = map_Ensembl_IDs_to_Entrez(snp_egene_el, ensembl_colname="gene_ensembl_id", new_colname="gene_id", keep_unmapped=True)
    snp_egene_el.sort_values(by="gene_id", inplace=True)

    output_dir = get_path("edge/snp-gene")
    snp_egene_el.to_csv(os.path.join(output_dir, "SNP_x_GTEx.edgelist"), sep="\t", index=False, header=False)
