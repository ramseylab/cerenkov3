import pandas as pd


_dei_set = set(pd.read_csv("./Entrez_discontinued_id.tsv", sep="\t", 
                           usecols=["Discontinued_GeneID"], 
                           dtype={"Discontinued_GeneID": int}).Discontinued_GeneID)

def filter_dei(df, id_col):
    global _dei_set

    found_dei = set(df[id_col]).intersection(_dei_set)

    if found_dei:
        print("Found Disctinued Entrez ID: {}".format(found_dei))

    flag_dei = df[id_col].isin(found_dei)

    return df.loc[~flag_dei, :]