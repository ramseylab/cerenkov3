import os
import pandas as pd
from util_path import get_path

_gene_dir = get_path("vertex/gene")
_ensembl_entrez_map = pd.read_csv(os.path.join(_gene_dir, "Ensembl_x_Entrez.tsv"), 
                                  sep="\t", dtype={"Ensembl_Gene_ID": str, "Entrez_Gene_ID": str})
# rename columns in order to lower the possibility of naming conflicts
_ensembl_entrez_map.rename(columns={"Ensembl_Gene_ID": "_Ensembl_Gene_ID", 
                                    "Entrez_Gene_ID": "_Entrez_Gene_ID"}, inplace=True)

def get_mapped_Ensembl_IDs():
    return set(_ensembl_entrez_map._Ensembl_Gene_ID)

def map_Ensembl_IDs_to_Entrez(df, ensembl_colname, new_colname, keep_unmapped=False):
    """
    Please make sure that no column in `df` is named "_Ensembl_Gene_ID" or "_Entrez_Gene_ID"
    """
    mapped_df = df.merge(_ensembl_entrez_map, how="left", left_on=ensembl_colname, right_on="_Ensembl_Gene_ID")
    mapped_df.drop("_Ensembl_Gene_ID", axis=1, inplace=True)
    mapped_df.rename(columns={ensembl_colname: "_Ensembl_Gene_ID"}, inplace=True)

    if keep_unmapped:  # keep unmapped Ensembl IDs
        def first_non_nan(x):
            if not pd.isnull(x["_Entrez_Gene_ID"]):
                return x["_Entrez_Gene_ID"]
            
            return x["_Ensembl_Gene_ID"]

        rename_kwargs = {new_colname: mapped_df.apply(first_non_nan, axis=1)}
        mapped_df = mapped_df.assign(**rename_kwargs)
        mapped_df.drop(["_Ensembl_Gene_ID", "_Entrez_Gene_ID"], axis=1, inplace=True)

        return mapped_df.drop_duplicates()
    else:
        mapped_df.dropna(inplace=True)
        mapped_df.rename(columns={"_Entrez_Gene_ID": new_colname}, inplace=True)
        mapped_df.drop("_Ensembl_Gene_ID", axis=1, inplace=True)

        return mapped_df.drop_duplicates()
