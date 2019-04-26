import os
import pandas as pd
from util_path import get_path

if __name__ == "__main__":
    hn_dir = get_path("edge/gene-gene/co-expression/HumanNet")
    humannet = pd.read_csv(os.path.join(hn_dir, "HumanNet-XN.tsv"), sep="\t",
                           comment="#", names=["Gene_A", "Gene_B", "LLS"], usecols=["Gene_A", "Gene_B"])

    humannet = humannet.drop_duplicates().sort_values(by=["Gene_A", "Gene_B"])
    humannet.to_csv(os.path.join(hn_dir, "p1_reduced_HumanNet.tsv"), sep="\t", index=False, header=True)
