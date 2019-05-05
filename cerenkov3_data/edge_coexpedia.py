import os
import pandas as pd
from functools import reduce
from util_path import get_path
from util_edge import remove_duplicated_undirected_edges


def read_coexpedia_dfms(_dir):
    for fn in os.listdir(_dir):
        path = os.path.join(_dir, fn)        
        dfm = pd.read_csv(path, sep="\t", header=None, names=["Gene_A", "Gene_B", "LLS"])
        
        geo_accession_id = fn.split(".txt")[0]
        dfm = dfm.assign(GEO_ACCESSION=geo_accession_id)

        yield dfm


if __name__ == "__main__":
    res_dir = get_path("resource/Coexpedia")
    dfms = read_coexpedia_dfms(os.path.join(res_dir, "Hsa"))

    merged_dfm = pd.concat(dfms, axis=0)
    merged_dfm.to_csv(os.path.join(res_dir, "p1_coexpedia_merged.tsv"), sep="\t", index=False)
    
    reduced_dfm = merged_dfm.loc[:, ["Gene_A", "Gene_B"]]
    reduced_dfm = remove_duplicated_undirected_edges(reduced_dfm, sort=True)
    
    edge_dir = get_path("edge/gene-gene")
    reduced_dfm.to_csv(os.path.join(edge_dir, "Coexpedia.edgelist"), sep="\t", index=False, header=False)