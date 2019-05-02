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
    input_dir = get_path("resource/Coexpedia/Hsa")
    dfms = read_coexpedia_dfms(input_dir)
    
    output_dir = get_path("edge/gene-gene/co-expression")

    merged_dfm = pd.concat(dfms, axis=0)
    merged_dfm.to_csv(os.path.join(output_dir, "coexpedia_p1_merged.tsv"), sep="\t", index=False)

    reduced_dfm = merged_dfm.loc[:, ["Gene_A", "Gene_B"]].drop_duplicates().sort_values(by=["Gene_A", "Gene_B"])
    reduced_dfm.to_csv(os.path.join(output_dir, "coexpedia.edgelist"), sep="\t", index=False, header=False)