#!/bin/bash

wget 'https://4dgenome.research.chop.edu/Tables/4DGenome_HomoSapiens_hg19.txt'
mv 4DGenome_HomoSapiens_hg19.txt ./edge/snp-gene/Chromatin-Interaction/4DGenome/

Rscript edge_4DGenome_preprocess.R
Rscript edge_4DGenome_make_BED.R

# Run bedtools intersect

## Modify this path if you have your own SNP bed file 
snp_bed="./vertex/SNP/OSU18/osu18_SNP.bed"

fd_genome_dir="./edge/snp-gene/Chromatin-Interaction/4DGenome"
ensembl_dir="./vertex/gene/Ensembl"

bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorA.bed -b ${snp_bed} -wo > ${fd_genome_dir}/4DGenome_InteractorA_SNP.bed
bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorB.bed -b ${snp_bed} -wo > ${fd_genome_dir}/4DGenome_InteractorB_SNP.bed

bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorA.bed -b ${ensembl_dir}/p2_Ensembl_TSS.bed -wo > ${fd_genome_dir}/4DGenome_InteractorA_Ensembl_TSS.bed
bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorB.bed -b ${ensembl_dir}/p2_Ensembl_TSS.bed -wo > ${fd_genome_dir}/4DGenome_InteractorB_Ensembl_TSS.bed

bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorA.bed -b ${ensembl_dir}/p2_Ensembl_promoter.bed -wo > ${fd_genome_dir}/4DGenome_InteractorA_Ensembl_promoter.bed
bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorB.bed -b ${ensembl_dir}/p2_Ensembl_promoter.bed -wo > ${fd_genome_dir}/4DGenome_InteractorB_Ensembl_promoter.bed

python3 edge_4DGenome.py