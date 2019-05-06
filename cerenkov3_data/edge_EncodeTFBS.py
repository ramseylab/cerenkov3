import os
import pandas as pd
import numpy as np
import mygene
from util_uscs import GenomeBrowserClient as GBC
from util_uscs import binary_encode_tfbs
from util_path import get_path
from util_dei import filter_dei


def get_snp_tfbs_matrix(rsid_list, config_key="remote_hg19"):
    # config_key == "remote_hg19" or "local_hg19"
    with GBC(config_key) as gbc:
        _temp_df = gbc.select_tfbs(rsid_list)
        snp_tfbs_matrix = binary_encode_tfbs(_temp_df, target_colname="tfName", value_sep=',', dest_colname_prefix=None).set_index("name") 

    return snp_tfbs_matrix

def update_tfbs_symbol(snp_tfbs_matrix):
    # ----- Get rid of deprecated/informal symbols ----- #

    # --- Rename --- #
        # FAM48A => SUPT20H; see https://www.ncbi.nlm.nih.gov/gene/55578
        # RDBP => NELFE; see https://www.ncbi.nlm.nih.gov/gene/7936
        # RPC155 => POLR3A; see https://www.ncbi.nlm.nih.gov/gene/11128
        # SREBP1 => SREBF1; see https://www.ncbi.nlm.nih.gov/gene/6720

    # --- Collapse --- #
        # KAP1 => TRIM28; see https://www.ncbi.nlm.nih.gov/gene/10155
        # GRp20 => GR => NR3C1; see https://www.ncbi.nlm.nih.gov/gene/2908
        #   - GRp20 == GR (P-20)
        # SIN3AK20 => SIN3A; see https://www.ncbi.nlm.nih.gov/gene/25942
        #   - SIN3AK20 == SIN3A (K-20)
    
    if "FAM48A" in snp_tfbs_matrix.columns:
        snp_tfbs_matrix.rename(columns={"FAM48A": "SUPT20H"}, inplace=True)
    if "RDBP" in snp_tfbs_matrix.columns:
        snp_tfbs_matrix.rename(columns={"RDBP": "NELFE"}, inplace=True)
    if "RPC155" in snp_tfbs_matrix.columns:
        snp_tfbs_matrix.rename(columns={"RPC155": "POLR3A"}, inplace=True)
    if "SREBP1" in snp_tfbs_matrix.columns:
        snp_tfbs_matrix.rename(columns={"SREBP1": "SREBF1"}, inplace=True)

    if "KAP1" in snp_tfbs_matrix.columns:
        if "TRIM28" in snp_tfbs_matrix.columns:
            snp_tfbs_matrix.loc[:, "TRIM28"] |= snp_tfbs_matrix.loc[:, "KAP1"] 
            snp_tfbs_matrix.drop(['KAP1'], axis=1, inplace=True)
        else:
            snp_tfbs_matrix.rename(columns={"KAP1": "TRIM28"}, inplace=True)
    if "GRp20" in snp_tfbs_matrix.columns:
        if "NR3C1" in snp_tfbs_matrix.columns:
            snp_tfbs_matrix.loc[:, "NR3C1"] |= snp_tfbs_matrix.loc[:, "GRp20"] 
            snp_tfbs_matrix.drop(['GRp20'], axis=1, inplace=True)
        else:
            snp_tfbs_matrix.rename(columns={"GRp20": "NR3C1"}, inplace=True)
    if "SIN3AK20" in snp_tfbs_matrix.columns:
        if "SIN3A" in snp_tfbs_matrix.columns:
            snp_tfbs_matrix.loc[:, "SIN3A"] |= snp_tfbs_matrix.loc[:, "SIN3AK20"] 
            snp_tfbs_matrix.drop(['SIN3AK20'], axis=1, inplace=True)
        else:
            snp_tfbs_matrix.rename(columns={"SIN3AK20": "SIN3A"}, inplace=True)
    
    # Restore the column names
    # snp_tfbs_matrix.rename(columns=lambda x: "tf_{}".format(x), inplace=True)

    return snp_tfbs_matrix

def convert_to_map(snp_tfbs_matrix, map_colnames):
    """
    Convert a binary matrix of SNPs by TFBSs to a long 2-column dataframe. E.g.  
    
    Input:
    
        name        ATF3    BACH1   BCL3    BHLHE40 CCNT2   CHD2
        rs1536167	   1        1	   0          1	    1	   1

    Output:

        name        tfbs
        rs1536167   ATF3
        rs1536167   BACH1
        rs1536167   BHLHE40
        rs1536167   ATFCCNT23
        rs1536167   CHD2

    The input binary matrix MUST use RSIDs as the only index and TFBS symbols as the only column names.
    """
    # Convert the binary matrix into a list of (rsid, TFBS) tuples
    snp_tfbs_map = [(snp_tfbs_matrix.index[row], snp_tfbs_matrix.columns[col]) for row, col in zip(*np.where(snp_tfbs_matrix==1))]
    # Convert the the list of tuples into a dataframe
    snp_tfbs_map = pd.DataFrame(snp_tfbs_map, columns=map_colnames)

    return snp_tfbs_map

def get_rsid_list():
    snp_dir = get_path("vertex/SNP")
    snp_fn = "osu18_SNP.bed"
    rsid_list = pd.read_csv(os.path.join(snp_dir, snp_fn), sep="\t", header=None, usecols=[3], 
                            names=["name"], squeeze=True).tolist()
    
    return rsid_list


if __name__ == "__main__":
    res_dir = get_path("resource/EncodeTFBS")

    rsid_list = get_rsid_list()
    snp_tfbs_matrix = get_snp_tfbs_matrix(rsid_list, config_key="local_hg19")
    snp_tfbs_matrix = update_tfbs_symbol(snp_tfbs_matrix)
    snp_tfbs_matrix.to_csv(os.path.join(res_dir, "p1_SNP_x_TFBS_matrix.tsv"), sep="\t", header=True, index=True)
    
    snp_tfbs_map = convert_to_map(snp_tfbs_matrix, map_colnames=["name", "symbol"])
    snp_tfbs_map.to_csv(os.path.join(res_dir, "p2_SNP_x_TFBS_map.tsv"), sep="\t", header=True, index=False)

    mg_api = mygene.MyGeneInfo()
    tfbs_list = snp_tfbs_matrix.columns.tolist()
    tfbs_gene_map = mg_api.querymany(tfbs_list, scopes='symbol', species=9606, as_dataframe=True)  # taxid 9606 => Homo sapiens 
    tfbs_gene_map = tfbs_gene_map.loc[:, ["entrezgene"]]  # TFBSs are stored in the index of this dataframe
    tfbs_gene_map.index.name = 'symbol'
    tfbs_gene_map = tfbs_gene_map.reset_index().sort_values(by="symbol")
    tfbs_gene_map.to_csv(os.path.join(res_dir, "p3_TFBS_x_Gene_map.tsv"), sep="\t", header=True, index=False)

    snp_gene_map = snp_tfbs_map.merge(tfbs_gene_map, on="symbol")
    snp_gene_el = snp_gene_map.loc[:, ["name", "entrezgene"]]
    snp_gene_el = filter_dei(snp_gene_el, id_col="entrezgene")

    output_dir = get_path("edge/snp-gene/")
    snp_gene_el.to_csv(os.path.join(output_dir, "SNP_x_EncodeTFBS.edgelist"), sep="\t", header=False, index=False)