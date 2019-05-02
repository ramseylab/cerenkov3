import os
import pandas as pd
from util_path import get_path

if __name__ == "__main__":
    input_dir = get_path("resource/HumanNet")
    humannet = pd.read_csv(os.path.join(input_dir, "HumanNet-XN.tsv"), sep="\t",
                           comment="#", names=["Gene_A", "Gene_B", "LLS"], usecols=["Gene_A", "Gene_B"])

    output_dir = get_path("edge/gene-gene/co-expression")
    humannet = humannet.drop_duplicates().sort_values(by=["Gene_A", "Gene_B"])
    humannet.to_csv(os.path.join(output_dir, "HumanNet.edgelist"), sep="\t", index=False, header=False)
