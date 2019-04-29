#!/bin/bash

wget 'https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz'

tar -xzf GTEx_Analysis_v7_eQTL.tar.gz
rm GTEx_Analysis_v7_eQTL.tar.gz

gunzip GTEx_Analysis_v7_eQTL/*.egenes.txt.gz
rm GTEx_Analysis_v7_eQTL/*.signif_variant_gene_pairs.txt.gz
mv GTEx_Analysis_v7_eQTL ./edge/snp-gene/eQTL/GTEx/

python3 edge_GTEx.py