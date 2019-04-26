#!/bin/bash

wget 'https://4dgenome.research.chop.edu/Tables/4DGenome_HomoSapiens_hg19.txt'
mv 4DGenome_HomoSapiens_hg19.txt ./edge/snp-gene/Chromatin-Interaction/4DGenome/

Rscript edge_4DGenome_preprocess.R
Rscript edge_4DGenome_make_BED.R

fd_genome_dir="./edge/snp-gene/Chromatin-Interaction/4DGenome"
osu18_snp_dir="./vertex/SNP/OSU18"
ensembl_dir="./vertex/gene/Ensembl"

bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorA.bed -b ${osu18_snp_dir}/osu18_SNP.bed -wo > ${fd_genome_dir}/InteractorA_SNP.bed
bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorB.bed -b ${osu18_snp_dir}/osu18_SNP.bed -wo > ${fd_genome_dir}/InteractorB_SNP.bed

bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorA.bed -b ${ensembl_dir}/ensembl_gene_TSS.bed -wo > ${fd_genome_dir}/InteractorA_TSS.bed
bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorB.bed -b ${ensembl_dir}/ensembl_gene_TSS.bed -wo > ${fd_genome_dir}/InteractorB_TSS.bed

bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorA.bed -b ${ensembl_dir}/ensembl_gene_promoter.bed -wo > ${fd_genome_dir}/InteractorA_promoter.bed
bedtools intersect -a ${fd_genome_dir}/4DGenome_InteractorB.bed -b ${ensembl_dir}/ensembl_gene_promoter.bed -wo > ${fd_genome_dir}/InteractorB_promoter.bed
