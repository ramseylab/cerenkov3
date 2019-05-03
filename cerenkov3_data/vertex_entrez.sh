#!/bin/bash

DIRECTORY="./resource/Entrez"
G2E_FILENAME="gene2ensembl_9606.tsv"
BIOMART_G37_FILENAME="GRCh37_p13_mart_export.txt"
BIOMART_G38_FILENAME="GRCh38_p12_mart_export.txt"
GENCODE_G37_FILENAME="gencode.v29lift37.metadata.EntrezGene"
GENCODE_G38_FILENAME="gencode.v29.metadata.EntrezGene"
HGNC_FILENAME="HGNC_custom.txt"

# gene2ensembl

if [ -f "${DIRECTORY}/${G2E_FILENAME}" ]; then
    echo "[Entrez] file '${DIRECTORY}/${G2E_FILENAME}' exists; no downloading"
else
    wget 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz'
    # Extract rows of Human genes only (whose Ensembl ID starts with "ENSG0")
    zcat gene2ensembl.gz | grep ENSG0 | cut -f 1,2,3 > ${G2E_FILENAME}
    rm gene2ensembl.gz
    mv ${G2E_FILENAME} ${DIRECTORY}
fi

# BioMart

if [ -f "${DIRECTORY}/${BIOMART_G37_FILENAME}" ]; then
    echo "[Entrez] file '${DIRECTORY}/${BIOMART_G37_FILENAME}' exists; no downloading"
else
    wget -O ${BIOMART_G37_FILENAME} 'http://grch37.ensembl.org/biomart/martservice/results?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query> <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="ensembl_transcript_id"/><Attribute name="hgnc_id"/></Dataset></Query>'
    mv ${BIOMART_G37_FILENAME} ${DIRECTORY}
fi

if [ -f "${DIRECTORY}/${BIOMART_G38_FILENAME}" ]; then
    echo "[Entrez] file '${DIRECTORY}/${BIOMART_G38_FILENAME}' exists; no downloading"
else
    wget -O ${BIOMART_G38_FILENAME} 'http://apr2019.archive.ensembl.org/biomart/martservice/results?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id" /><Attribute name="ensembl_transcript_id" /><Attribute name="hgnc_id" /></Dataset></Query>'
    mv ${BIOMART_G38_FILENAME} ${DIRECTORY}
fi

# Gencode

if [ -f "${DIRECTORY}/${GENCODE_G37_FILENAME}" ]; then
    echo "[Entrez] file '${DIRECTORY}/${GENCODE_G37_FILENAME}' exists; no downloading"
else
    wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/${GENCODE_G37_FILENAME}.gz"
    gunzip ${GENCODE_G37_FILENAME}.gz
    mv ${GENCODE_G37_FILENAME} ${DIRECTORY}
fi

if [ -f "${DIRECTORY}/${GENCODE_G38_FILENAME}" ]; then
    echo "[Entrez] file '${DIRECTORY}/${GENCODE_G38_FILENAME}' exists; no downloading"
else
    wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/${GENCODE_G38_FILENAME}.gz"
    gunzip ${GENCODE_G38_FILENAME}.gz
    mv ${GENCODE_G38_FILENAME} ${DIRECTORY}
fi

# HGNC

if [ -f "${DIRECTORY}/${HGNC_FILENAME}" ]; then
    echo "[Entrez] file '${DIRECTORY}/${HGNC_FILENAME}' exists; no downloading"
else
    wget -O ${HGNC_FILENAME} 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_chrom_map&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=family.id&col=family.name&col=md_eg_id&col=md_ensembl_id&col=md_refseq_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'
    mv ${HGNC_FILENAME} ${DIRECTORY}
fi

# Map Ensembl IDs to Entrez; convert the IDs in Ensembl BED files to Entrez 

python3 vertex_entrez.py

echo "[Entrez] done!"