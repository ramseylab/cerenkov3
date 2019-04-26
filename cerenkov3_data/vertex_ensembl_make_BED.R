# This script aims to:
    # 1. use `p1_ensembl_gene_df.tsv` as input
    # 2. save "chromosome", "promoter_start", "promoter_end", "gene_name", "ensembl_gene_id", "entrez_gene_id" and "strand" to a BED file, `ensembl_gene_promoter.bed`
    #     - You always have `transcript_start < transcript_end`, `gene_start < gene_end` and `TSS_start + 1 == TSS_end`
    #     - If strand == +1, define:
    #         - promoter_start = TSS_start - 2000, and 
    #         - promoter_end = TSS_end + 500 (smaller coordinates mean upstream)
    #     - If strand == -1, define:
    #         - promoter_start = TSS_start - 500, and 
    #         - promoter_end = TSS_end + 2000 (bigger coordinates mean upstream)
    # 3. save "chromosome", "TSS_start", "TSS_end", "gene_name", "ensembl_gene_id", "entrez_gene_id" and "strand" to a BED file, `ensembl_gene_TSS.bed`
    # 4. Note: BED format is not strictly followed here.

library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

ensembl_dir <- "./vertex/gene/Ensembl/"

ensembl_df <- read_tsv(paste(ensembl_dir, "p1_ensembl_gene_df.tsv", sep=""))

ensembl_df %>% 
	mutate(promoter_start = case_when(strand == 1 ~ TSS_start - 2000, 
                     	              strand == -1 ~ TSS_start - 500), 
		   promoter_end = case_when(strand == 1 ~ TSS_end + 500, 
                     	              strand == -1 ~ TSS_end + 2000)) %>%
	select(chromosome, promoter_start, promoter_end, gene_name, ensembl_gene_id, entrez_gene_id, strand) %>%
	mutate(entrez_gene_id = replace_na(entrez_gene_id, ".")) %>%
	distinct(.keep_all = TRUE) %>% 
	write.table(paste(ensembl_dir, "p2_ensembl_gene_promoter.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

ensembl_df %>% 
	select(chromosome, TSS_start, TSS_end, gene_name, ensembl_gene_id, entrez_gene_id, strand) %>%
	mutate(entrez_gene_id = replace_na(entrez_gene_id, ".")) %>%
	distinct(.keep_all = TRUE) %>% 
	write.table(paste(ensembl_dir, "p2_ensembl_gene_TSS.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

## save "chromosome", "transcript_start", "transcript_end", "gene_name", "ensembl_gene_id", "entrez_gene_id" and "strand" to a BED file, `ensembl_gene_transcript.bed`
# ensembl_df %>% 
# 	select(chromosome, transcript_start, transcript_end, gene_name, ensembl_gene_id, entrez_gene_id, strand) %>%
# 	mutate(entrez_gene_id = replace_na(entrez_gene_id, ".")) %>% 
# 	distinct(.keep_all = TRUE) %>% 
# 	write.table("p2_ensembl_gene_transcript.bed", sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)

## save "chromosome", "gene_start", "gene_end", "gene_name", "ensembl_gene_id", "entrez_gene_id" and "strand" to a BED file, `ensembl_gene_itself.bed`
# ensembl_df %>% 
# 	select(chromosome, gene_start, gene_end, gene_name, ensembl_gene_id, entrez_gene_id, strand) %>%
# 	mutate(entrez_gene_id = replace_na(entrez_gene_id, ".")) %>%
# 	distinct(.keep_all = TRUE) %>% 
# 	write.table("p2_ensembl_gene_itself.bed", sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)