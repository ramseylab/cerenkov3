import os
import pandas as pd
from util_path import get_path

_dei_dir = get_path("util/discontinued_entrez_id")
_dei_set = set(pd.read_csv(os.path.join(_dei_dir, "discontinued_entrez_id.tsv"), 
                           sep="\t", usecols=["Discontinued_GeneID"], 
                           dtype={"Discontinued_GeneID": int}).Discontinued_GeneID)

def filter_dei(df, id_col):
    global _dei_set

    found_dei = set(df[id_col]).intersection(_dei_set)

    if found_dei:
        print("Found Disctinued Entrez ID: {}".format(found_dei))

    flag_dei = df[id_col].isin(found_dei)

    return df.loc[~flag_dei, :]