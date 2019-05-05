#!/bin/bash

DIRECTORY="./resource/GTEx"
SUB_DIR="GTEx_Analysis_v7_eQTL"

if [ -d "${DIRECTORY}/${SUB_DIR}" ]; then
    # Control will enter here if that folder exists.
    echo "[GTEx] folder '${DIRECTORY}/${SUB_DIR}' exists; skip downloading"
else
    wget "https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/${SUB_DIR}.tar.gz"
    tar -xzf ${SUB_DIR}.tar.gz
    rm ${SUB_DIR}.tar.gz

    gunzip ${SUB_DIR}/*.egenes.txt.gz
    rm ${SUB_DIR}/*.signif_variant_gene_pairs.txt.gz
    mv ${SUB_DIR} ${DIRECTORY}
fi

python3 edge_GTEx.py

echo "[GTEx] done!"