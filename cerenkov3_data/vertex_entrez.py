import os
import pandas as pd
import mygene
from util_path import get_path
from util_dei import filter_dei

res_dir = get_path("resource/Entrez")
gene_dir = get_path("vertex/gene")

mg = mygene.MyGeneInfo()

def read_gene2ensembl():
    global res_dir

    g2e_df = pd.read_csv(os.path.join(res_dir, "gene2ensembl_9606.tsv"), sep="\t", header=None, 
                         names=["Tax_ID", "Entrez_Gene_ID", "Ensembl_Gene_ID"])

    unique_tax_ids = g2e_df.Tax_ID.unique()
    assert(len(unique_tax_ids) == 1)
    assert(unique_tax_ids[0] == 9606)

    g2e_df = g2e_df.drop("Tax_ID", axis=1).drop_duplicates().reindex()

    return g2e_df

def read_biomart(GRCh):
    global res_dir

    if GRCh == "37":
        filename = "GRCh37_p13_mart_export.txt"
    elif GRCh == "38":
        filename = "GRCh38_p12_mart_export.txt"
    else:
        raise ValueError("Cannot recognize GRCh param. Got {}. Use string '37' or '38' instead.".format(GRCh))

    biomart_df = pd.read_csv(os.path.join(res_dir, filename), sep="\t", dtype=str)
    biomart_df = biomart_df.rename(columns={"Gene stable ID": "Ensembl_Gene_ID",
                                            "HGNC ID": "HGNC_ID", 
                                            "Transcript stable ID": "Ensembl_Tx_ID"})

    if GRCh == "37":
        return biomart_df
    else:  # GRCh == "38"
        # Note that in GRCh38 BioMart data, every HGNC ID starts with a prefix "HGNC:"
        # No such prefix in GRCh37 BioMart data
        if all(biomart_df.HGNC_ID.dropna().str.startswith("HGNC:")):
            biomart_df.loc[:, "HGNC_ID"] = biomart_df.loc[:, "HGNC_ID"].apply(lambda x: x.split(":")[1] if not pd.isnull(x) else x)
        
        return biomart_df

def read_gencode(GRCh):
    global res_dir

    if GRCh == "37":
        filename = "gencode.v29lift37.metadata.EntrezGene"
    elif GRCh == "38":
        filename = "gencode.v29.metadata.EntrezGene"
    else:
        raise ValueError("Cannot recognize GRCh param. Got {}. Use string '37' or '38' instead.".format(GRCh))

    gencode_df = pd.read_csv(os.path.join(res_dir, filename), sep="\t", header=None, 
                             names=["Ensembl_Tx_ID", "Entrez_Gene_ID_gc"], dtype=str)
    gencode_df.loc[:, "Ensembl_Tx_ID"] = gencode_df.loc[:, "Ensembl_Tx_ID"].apply(lambda x: x.split(".")[0])

    return gencode_df

def read_hgnc():
    global res_dir 

    hgnc_df = pd.read_csv(os.path.join(res_dir, "HGNC_custom.txt"), sep="\t", 
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

def _map_biomart_gencode_hgnc(biomart_df, gencode_df, hgnc_df):
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

def map_biomart_gencode_hgnc(GRCh):
    biomart_df = read_biomart(GRCh)
    gencode_df = read_gencode(GRCh)
    hgnc_df = read_hgnc()

    return _map_biomart_gencode_hgnc(biomart_df, gencode_df, hgnc_df)

def read_ensembl_bed(filename):
    global gene_dir

    bed_df = pd.read_csv(os.path.join(gene_dir, filename), sep="\t", usecols=[0,1,2,3,4,5], header=None, 
                         names=["chrom", "chromStart", "chromEnd", "name", "id", "strand"])

    return bed_df

def convert_ensembl_bed_to_entrez(ensembl_bed, ensembl_entrez_map):
    merged_df = ensembl_entrez_map.merge(ensembl_bed, how="inner", left_on="Ensembl_Gene_ID", right_on="id")
    entrez_bed = merged_df.loc[:, ["chrom", "chromStart", "chromEnd", "name", "Entrez_Gene_ID", "strand"]].drop_duplicates()
    entrez_bed.sort_values(by=["chrom", "chromStart"], inplace=True)

    return entrez_bed


if __name__ == "__main__":
    # Mapping Ensembl IDs to Entrez
    g2e_map = read_gene2ensembl()
    grch37_map = map_biomart_gencode_hgnc(GRCh="37")
    grch38_map = map_biomart_gencode_hgnc(GRCh="38")

    g2e_map.to_csv(os.path.join(res_dir, "p1_Ensembl_x_Entrez_gene2ensembl.tsv"), sep="\t", header=True, index=False)
    grch37_map.to_csv(os.path.join(res_dir, "p1_Ensembl_x_Entrez_grch37_biomart.tsv"), sep="\t", header=True, index=False)
    grch38_map.to_csv(os.path.join(res_dir, "p1_Ensembl_x_Entrez_grch38_biomart.tsv"), sep="\t", header=True, index=False)

    ensembl_entrez_map = pd.concat([grch37_map, grch38_map, g2e_map], axis=0, sort=True).drop_duplicates()
    ensembl_entrez_map = filter_dei(ensembl_entrez_map, "Entrez_Gene_ID")
    ensembl_entrez_map.to_csv(os.path.join(gene_dir, "Ensembl_x_Entrez.tsv"), sep="\t", header=True, index=False)

    # Convert the IDs in Ensembl BED files to Entrez
    ensembl_prm_bed = read_ensembl_bed("Ensembl_promoter.bed")
    ensembl_tss_bed = read_ensembl_bed("Ensembl_TSS.bed")

    entrez_prm_bed = convert_ensembl_bed_to_entrez(ensembl_prm_bed, ensembl_entrez_map)
    entrez_tss_bed = convert_ensembl_bed_to_entrez(ensembl_tss_bed, ensembl_entrez_map)

    entrez_prm_bed.to_csv(os.path.join(gene_dir, "Entrez_promoter.bed"), sep="\t", header=False, index=False)
    entrez_tss_bed.to_csv(os.path.join(gene_dir, "Entrez_TSS.bed"), sep="\t", header=False, index=False)