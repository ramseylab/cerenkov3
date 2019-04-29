# This script aims to:
    # 1. use `p1_ensembl_gene_df.tsv` as input
    # 2. save "chromosome", "promoter_start", "promoter_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file
    #     - You always have `transcript_start < transcript_end`, `gene_start < gene_end` and `TSS_start + 1 == TSS_end`
    #     - If strand == +1, define:
    #         - promoter_start = TSS_start - 2000, and 
    #         - promoter_end = TSS_end + 500 (smaller coordinates mean upstream)
    #     - If strand == -1, define:
    #         - promoter_start = TSS_start - 500, and 
    #         - promoter_end = TSS_end + 2000 (bigger coordinates mean upstream)
    # 3. save "chromosome", "TSS_start", "TSS_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file
# Note: BED format is not strictly followed here.

library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

ensembl_dir <- "./vertex/gene/Ensembl/"
# Note that this the coordinates in this file is already 0-based (converted by vertex_ensembl_preprocess.R) and compatible with BED format
ensembl_df <- read_tsv(paste(ensembl_dir, "p1_Ensembl.tsv", sep=""), 
					   col_types=cols(chromosome = col_character(),
								      transcript_start = col_integer(),
									  transcript_end = col_integer(),
									  strand = col_integer(),
									  ensembl_gene_id = col_character(),
									  gene_name = col_character(),
									  gene_start = col_integer(),
									  gene_end = col_integer(),
									  transcript_id = col_character(),
									  transcript_name = col_character(),
									  TSS_start = col_integer(),
									  TSS_end = col_integer()))

ensembl_df %>% 
	mutate(promoter_start = case_when(strand == 1 ~ TSS_start - 2000, 
                     	              strand == -1 ~ TSS_start - 500), 
		   promoter_end = case_when(strand == 1 ~ TSS_end + 500, 
                     	              strand == -1 ~ TSS_end + 2000)) %>%
	select(chromosome, promoter_start, promoter_end, gene_name, ensembl_gene_id, strand) %>%
	distinct(.keep_all = TRUE) %>% 
	write.table(paste(ensembl_dir, "p2_Ensembl_promoter.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

ensembl_df %>% 
	select(chromosome, TSS_start, TSS_end, gene_name, ensembl_gene_id, strand) %>%
	distinct(.keep_all = TRUE) %>% 
	write.table(paste(ensembl_dir, "p2_Ensembl_TSS.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

# 4. save "chromosome", "transcript_start", "transcript_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file
# ensembl_df %>% 
# 	select(chromosome, transcript_start, transcript_end, gene_name, ensembl_gene_id, strand) %>%
# 	distinct(.keep_all = TRUE) %>% 
# 	write.table(paste(ensembl_dir, "p2_Ensembl_transcript.bed", sep=""), sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)

## save "chromosome", "gene_start", "gene_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file, `ensembl_gene_itself.bed`
# ensembl_df %>% 
# 	select(chromosome, gene_start, gene_end, gene_name, ensembl_gene_id, strand) %>%
# 	distinct(.keep_all = TRUE) %>% 
# 	write.table("p2_Ensembl_gene.bed", sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)