import os
import pandas as pd
from util_path import get_path
from util_edge import remove_duplicated_undirected_edges

if __name__ == "__main__":
    input_dir = get_path("resource/BioGRID")
    edgelist_df = pd.read_csv(os.path.join(input_dir, "BIOGRID-Human.tab2.tsv"), sep="\t", usecols=["Entrez Gene Interactor A", "Entrez Gene Interactor B"])

    edgelist_df = remove_duplicated_undirected_edges(edgelist_df, sort=True)

    output_dir = get_path("edge/gene-gene")
    edgelist_df.to_csv(os.path.join(output_dir, "BioGRID.edgelist"), sep="\t", index=False, header=False)
