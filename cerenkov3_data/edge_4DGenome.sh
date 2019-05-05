#!/bin/bash

DIRECTORY="./resource/4DGenome"
FILENAME="4DGenome_HomoSapiens_hg19.txt"

if [ -f "${DIRECTORY}/${FILENAME}" ]; then
    echo "[4DGenome] file '${DIRECTORY}/${FILENAME}' exists; skip downloading"
else
    wget "https://4dgenome.research.chop.edu/Tables/${FILENAME}"
    mv ${FILENAME} ${DIRECTORY}
fi

Rscript edge_4DGenome_preprocess.R

# Run bedtools intersect

snp_dir="./vertex/SNP"
snp_fn="osu18_SNP.bed"  ## Modify this filename if you have your own SNP bed file 
gene_dir="./vertex/gene"
tss_fn="Ensembl_TSS.bed"
prm_fn="Ensembl_promoter.bed"

bedtools intersect -a ${DIRECTORY}/InteractorA.bed -b ${snp_dir}/${snp_fn} -wo > ${DIRECTORY}/InteractorA_SNP_intxn.bed
bedtools intersect -a ${DIRECTORY}/InteractorB.bed -b ${snp_dir}/${snp_fn} -wo > ${DIRECTORY}/InteractorB_SNP_intxn.bed

bedtools intersect -a ${DIRECTORY}/InteractorA.bed -b ${gene_dir}/${tss_fn} -wo > ${DIRECTORY}/InteractorA_Ensembl_TSS_intxn.bed
bedtools intersect -a ${DIRECTORY}/InteractorB.bed -b ${gene_dir}/${tss_fn} -wo > ${DIRECTORY}/InteractorB_Ensembl_TSS_intxn.bed

bedtools intersect -a ${DIRECTORY}/InteractorA.bed -b ${gene_dir}/${prm_fn} -wo > ${DIRECTORY}/InteractorA_Ensembl_promoter_intxn.bed
bedtools intersect -a ${DIRECTORY}/InteractorB.bed -b ${gene_dir}/${prm_fn} -wo > ${DIRECTORY}/InteractorB_Ensembl_promoter_intxn.bed

# Extract Edges

python3 edge_4DGenome.py

echo "[4DGenome] done!"