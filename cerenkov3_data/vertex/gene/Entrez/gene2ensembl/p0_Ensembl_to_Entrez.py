import pandas as pd
import sys
sys.path.append('../../../../util/discontinued_entrez_id')
from dei_util import filter_dei

def read_g2e_df():
    g2e_df = pd.read_csv("gene2ensembl_9606.tsv", sep="\t", header=None, 
                         names=["Tax_ID", "Entrez_Gene_ID", "Ensembl_Gene_ID"],
                         usecols=["Tax_ID", "Entrez_Gene_ID", "Ensembl_Gene_ID"])

    unique_tax_ids = g2e_df.Tax_ID.unique()
    assert(len(unique_tax_ids) == 1)
    assert(unique_tax_ids[0] == 9606)

    g2e_df = g2e_df.drop("Tax_ID", axis=1).drop_duplicates().reindex()

    return g2e_df


if __name__ == "__main__":
    g2e_df = read_g2e_df()

    g2e_df = filter_dei(g2e_df, "Entrez_Gene_ID")

    g2e_df.to_csv("Ensembl_x_Entrez.tsv", sep="\t", index=False)