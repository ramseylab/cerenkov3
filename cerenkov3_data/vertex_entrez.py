import os
import pandas as pd
import mygene
from util_path import get_path
from util_dei import filter_dei

entrez_dir = get_path("vertex/gene/Entrez")

mg = mygene.MyGeneInfo()

def gene2ensembl_map():
    global entrez_dir

    g2e_df = pd.read_csv(os.path.join(entrez_dir, "gene2ensembl_9606.tsv"), sep="\t", header=None, 
                         names=["Tax_ID", "Entrez_Gene_ID", "Ensembl_Gene_ID"])

    unique_tax_ids = g2e_df.Tax_ID.unique()
    assert(len(unique_tax_ids) == 1)
    assert(unique_tax_ids[0] == 9606)

    g2e_df = g2e_df.drop("Tax_ID", axis=1).drop_duplicates().reindex()

    return g2e_df

def read_grch37_biomart_df():
    global entrez_dir

    biomart_df = pd.read_csv(os.path.join(entrez_dir, "GRCh37_p13_mart_export.txt"), sep="\t", dtype=str)
    biomart_df = biomart_df.rename(columns={"Gene stable ID": "Ensembl_Gene_ID",
                                            "HGNC ID": "HGNC_ID", 
                                            "Transcript stable ID": "Ensembl_Tx_ID"})
    return biomart_df

def read_grch38_biomart_df():
    global entrez_dir

    biomart_df = pd.read_csv(os.path.join(entrez_dir, "GRCh38_p12_mart_export.txt"), sep="\t", dtype=str)
    biomart_df = biomart_df.rename(columns={"Gene stable ID": "Ensembl_Gene_ID",
                                            "HGNC ID": "HGNC_ID", 
                                            "Transcript stable ID": "Ensembl_Tx_ID"})

    # Note that in GRCh38 BioMart data, every HGNC ID starts with a prefix "HGNC:"
    # No such prefix in GRCh37 BioMart data
    if all(biomart_df.HGNC_ID.dropna().str.startswith("HGNC:")):
        biomart_df.loc[:, "HGNC_ID"] = biomart_df.loc[:, "HGNC_ID"].apply(lambda x: x.split(":")[1] if not pd.isnull(x) else x)
    
    return biomart_df

def read_gencode_df(fn):
    global entrez_dir

    gencode_df = pd.read_csv(os.path.join(entrez_dir, fn), sep="\t", header=None, 
                             names=["Ensembl_Tx_ID", "Entrez_Gene_ID_gc"], dtype=str)
    gencode_df.loc[:, "Ensembl_Tx_ID"] = gencode_df.loc[:, "Ensembl_Tx_ID"].apply(lambda x: x.split(".")[0])

    return gencode_df

def read_grch37_gencode_df():
    return read_gencode_df("gencode.v29lift37.metadata.EntrezGene")

def read_grch38_gencode_df():
    return read_gencode_df("gencode.v29.metadata.EntrezGene")

def read_hgnc_df():
    global entrez_dir 

    hgnc_df = pd.read_csv(os.path.join(entrez_dir, "HGNC_custom.txt"), sep="\t", 
                          usecols=["HGNC ID", "NCBI Gene ID", "Ensembl gene ID", 
                                   "NCBI Gene ID(supplied by NCBI)", "Ensembl ID(supplied by Ensembl)"], dtype=str)
    hgnc_df = hgnc_df.rename(columns={"HGNC ID": "HGNC_ID",
                                      "Ensembl gene ID": "Ensembl_Gene_ID", 
                                      "NCBI Gene ID": "Entrez_Gene_ID",
                                      "Ensembl ID(supplied by Ensembl)": "Ensembl_Gene_ID_external",
                                      "NCBI Gene ID(supplied by NCBI)": "Entrez_Gene_ID_external"})
    hgnc_df.loc[:, "Entrez_Gene_ID"] = hgnc_df.loc[:, "Entrez_Gene_ID"].fillna(hgnc_df.loc[:, "Entrez_Gene_ID_external"])
    hgnc_df.loc[:, "Ensembl_Gene_ID"] = hgnc_df.loc[:, "Ensembl_Gene_ID"].fillna(hgnc_df.loc[:, "Ensembl_Gene_ID_external"])
    hgnc_df = hgnc_df.drop(["Entrez_Gene_ID_external", "Ensembl_Gene_ID_external"], axis=1)

    hgnc_df = hgnc_df.rename(columns={"Entrez_Gene_ID": "Entrez_Gene_ID_hgnc"})

    hgnc_df.loc[:, "HGNC_ID"] = hgnc_df.loc[:, "HGNC_ID"].apply(lambda x: x.split(":")[1])

    return hgnc_df

def query_mygene(ensembl_ids):
    global mg

    mygene_df = mg.getgenes(ensembl_ids, field="entrezgene", species=9606, as_dataframe=True, verbose=False)

    mygene_df = mygene_df.loc[:, ["entrezgene"]]
    mygene_df = mygene_df.reset_index().rename(columns={"query": "Ensembl_Gene_ID", "entrezgene": "Entrez_Gene_ID_mg"})

    return mygene_df

def biomart_map(biomart_df, gencode_df, hgnc_df):
    merged_df = biomart_df.merge(gencode_df, how="left", on="Ensembl_Tx_ID")
    
    merged_df = merged_df.merge(hgnc_df.loc[:, ["Ensembl_Gene_ID", "Entrez_Gene_ID_hgnc"]], how="left", on="Ensembl_Gene_ID")
    merged_df = merged_df.merge(hgnc_df.loc[:, ["HGNC_ID", "Entrez_Gene_ID_hgnc"]], how="left", on="HGNC_ID", suffixes=["_1", "_2"])

    mygene_df = query_mygene(merged_df.Ensembl_Gene_ID.unique())
    merged_df = merged_df.merge(mygene_df, how="left", on="Ensembl_Gene_ID")

    ordered_keys = ["Entrez_Gene_ID_hgnc_2", "Entrez_Gene_ID_hgnc_1", "Entrez_Gene_ID_mg", "Entrez_Gene_ID_gc"]

    # drop a row if all its 4 "Entrez_Gene_ID_xx" columns are NA
    merged_df = merged_df.dropna(axis=0, how='all', subset=ordered_keys)
    merged_df = merged_df.drop(["Ensembl_Tx_ID", "HGNC_ID"], axis=1)
    merged_df = merged_df.drop_duplicates()

    def first_non_nan(row):
        nonlocal ordered_keys

        # next(iterator, default): Retrieve the next item from the iterator by calling its __next__() method. If default is given, it is returned if the iterator is exhausted, otherwise StopIteration is raised.
        return next((row[key] for key in ordered_keys if not pd.isnull(row[key])), None)
        
    merged_df = merged_df.assign(Entrez_Gene_ID = merged_df.apply(first_non_nan, axis=1).astype(int))

    merged_df = merged_df.drop(ordered_keys, axis=1)
    merged_df = merged_df.drop_duplicates()

    return merged_df

def grch37_biomart_map():
    biomart_df = read_grch37_biomart_df()
    gencode_df = read_grch37_gencode_df()
    hgnc_df = read_hgnc_df()

    return biomart_map(biomart_df, gencode_df, hgnc_df)

def grch38_biomart_map():
    biomart_df = read_grch38_biomart_df()
    gencode_df = read_grch38_gencode_df()
    hgnc_df = read_hgnc_df()

    return biomart_map(biomart_df, gencode_df, hgnc_df)

if __name__ == "__main__":
    g2e_df = gene2ensembl_map()
    grch37_df = grch37_biomart_map()
    grch38_df = grch38_biomart_map()

    g2e_df.to_csv(os.path.join(entrez_dir, "p1_Ensembl_x_Entrez_gene2ensembl.tsv"), sep="\t", header=True, index=False)
    grch37_df.to_csv(os.path.join(entrez_dir, "p1_Ensembl_x_Entrez_grch37_biomart.tsv"), sep="\t", header=True, index=False)
    grch38_df.to_csv(os.path.join(entrez_dir, "p1_Ensembl_x_Entrez_grch38_biomart.tsv"), sep="\t", header=True, index=False)

    merged_df = pd.concat([grch37_df, grch38_df, g2e_df], axis=0, sort=True).drop_duplicates()
    merged_df = filter_dei(merged_df, "Entrez_Gene_ID")

    merged_df.to_csv(os.path.join(entrez_dir, "p2_Ensembl_x_Entrez.tsv"), sep="\t", header=True, index=False)