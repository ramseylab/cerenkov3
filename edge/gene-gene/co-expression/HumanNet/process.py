import pandas as pd

if __name__ == "__main__":
    humannet = pd.read_csv("HumanNet-XN.tsv", sep="\t",
                           comment="#", names=["Gene_A", "Gene_B", "LLS"], usecols=["Gene_A", "Gene_B"])

    humannet = humannet.drop_duplicates().sort_values(by=["Gene_A", "Gene_B"])
    humannet.to_csv("p1_reduced_HumanNet.tsv", sep="\t", index=False, header=True)
