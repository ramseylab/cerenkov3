import pandas as pd
import os
from functools import reduce
from util_path import get_path


def read_coexpedia_dfms(_dir):
    for fn in os.listdir(_dir):
        path = os.path.join(_dir, fn)        
        dfm = pd.read_csv(path, sep="\t", header=None, names=["Gene_A", "Gene_B", "LLS"])
        
        geo_accession_id = fn.split(".txt")[0]
        dfm = dfm.assign(GEO_ACCESSION=geo_accession_id)

        yield dfm


if __name__ == "__main__":
    hsa_dir = get_path("edge/gene-gene/co-expression/Coexpedia/Hsa")
    dfms = read_coexpedia_dfms(hsa_dir)
    
    coexp_dir = get_path("edge/gene-gene/co-expression/Coexpedia")

    merged_dfm = pd.concat(dfms, axis=0)
    merged_dfm.to_csv(os.path.join(coexp_dir, "p1_merged_coexpedia.tsv"), sep="\t", index=False)

    reduced_dfm = merged_dfm.loc[:, ["Gene_A", "Gene_B"]].drop_duplicates().sort_values(by=["Gene_A", "Gene_B"])
    reduced_dfm.to_csv(os.path.join(coexp_dir, "p2_reduced_coexpedia.tsv"), sep="\t", index=False)