# This script aims to:
# 	1. read `mart_export.txt` and rename its columns
# 	2. since Ensembl uses 1-based coordinate system while we use 0-based, revise `TSS`, `gene_start` and `transcript_start` columns
#	3. add "chr" prefix to "chromosome" column, e.g. "1" => "chr1"
# 	4. save the revised dataframe to a TSV file
#   5. save "chromosome", "promoter_start", "promoter_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file
#      - You always have `transcript_start < transcript_end`, `gene_start < gene_end` and `TSS_start + 1 == TSS_end`
#      - If strand == +1, define:
#          - promoter_start = TSS_start - 2000, and 
#          - promoter_end = TSS_end + 500 (smaller coordinates mean upstream)
#      - If strand == -1, define:
#          - promoter_start = TSS_start - 500, and 
#          - promoter_end = TSS_end + 2000 (bigger coordinates mean upstream)
#   6. save "chromosome", "TSS_start", "TSS_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file
# Note: BED format is not strictly followed here.

suppressPackageStartupMessages(library(dplyr))
library(magrittr)
library(readr)

resource_dir <- "./resource/Ensembl/"

ensembl_gene_df <- read_tsv(paste(resource_dir, "mart_export.txt", sep=""), col_types=cols(
	`Chromosome/scaffold name` = col_character(),
	`Transcript start (bp)` = col_integer(),
	`Transcript end (bp)` = col_integer(),
	`Transcription start site (TSS)` = col_integer(),
	Strand = col_integer(),
	`Gene stable ID` = col_character(),
	`Gene name` = col_character(),
	`Gene start (bp)` = col_integer(),
	`Gene end (bp)` = col_integer(),
	`Transcript stable ID` = col_character(),
	`Transcript name` = col_character()
))

# Rename columns
ensembl_gene_df <- ensembl_gene_df %>% rename(
	ensembl_gene_id = `Gene stable ID`,
	gene_name = `Gene name`,
	strand = Strand,
	chromosome = `Chromosome/scaffold name`,
	TSS = `Transcription start site (TSS)`,
	transcript_id = `Transcript stable ID`, 
	transcript_name = `Transcript name`, 
	transcript_start = `Transcript start (bp)`,
	transcript_end = `Transcript end (bp)`,
	gene_start = `Gene start (bp)`,
	gene_end = `Gene end (bp)`
)

# Check if exceeded a integer's range
	# Otherwise you should use `col_double()` when reading the TSV
stopifnot(ensembl_gene_df$transcript_start %>% max < .Machine$integer.max)
stopifnot(ensembl_gene_df$transcript_end %>% max < .Machine$integer.max)
stopifnot(ensembl_gene_df$gene_start %>% max < .Machine$integer.max)
stopifnot(ensembl_gene_df$gene_end %>% max < .Machine$integer.max)

# "1" => "chr1", etc.
ensembl_gene_df <- ensembl_gene_df %>% mutate(chromosome = paste0("chr", chromosome))

# 1-based => 0-based
# See http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
# > To convert to 0-start, half-open: subtract 1 from start, end = same
ensembl_gene_df <- ensembl_gene_df %>% mutate(transcript_start = transcript_start - 1, 
						   					  gene_start = gene_start - 1)
# Extra to TSS: separate column "TSS" into two, e.g. TSS == 100 => TSS_start = 99; TSS_end = 100
ensembl_gene_df <- ensembl_gene_df %>% mutate(TSS_start = TSS - 1, 
											  TSS_end = TSS)

ensembl_gene_df %>% 
	select(-TSS) %>%  # Delete original `TSS` column				
	write_tsv(paste(resource_dir, "p1_Ensembl.tsv", sep=""))  # Save the revised dataframe

vertex_dir <- "./vertex/gene/"

ensembl_gene_df %>% 
	mutate(promoter_start = case_when(strand == 1 ~ TSS_start - 2000, 
                     	              strand == -1 ~ TSS_start - 500), 
		   promoter_end = case_when(strand == 1 ~ TSS_end + 500, 
                     	              strand == -1 ~ TSS_end + 2000)) %>%
	select(chromosome, promoter_start, promoter_end, gene_name, ensembl_gene_id, strand) %>%
	distinct(.keep_all = TRUE) %>% 
	write.table(paste(vertex_dir, "Ensembl_promoter.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

ensembl_gene_df %>% 
	select(chromosome, TSS_start, TSS_end, gene_name, ensembl_gene_id, strand) %>%
	distinct(.keep_all = TRUE) %>% 
	write.table(paste(vertex_dir, "Ensembl_TSS.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

# 7. save "chromosome", "transcript_start", "transcript_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file
# ensembl_gene_df %>% 
# 	select(chromosome, transcript_start, transcript_end, gene_name, ensembl_gene_id, strand) %>%
# 	distinct(.keep_all = TRUE) %>% 
# 	write.table(paste(vertex_dir, "Ensembl_transcript.bed", sep=""), sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)

# 8. save "chromosome", "gene_start", "gene_end", "gene_name", "ensembl_gene_id" and "strand" to a BED file, `ensembl_gene_itself.bed`
# ensembl_gene_df %>% 
# 	select(chromosome, gene_start, gene_end, gene_name, ensembl_gene_id, strand) %>%
# 	distinct(.keep_all = TRUE) %>% 
# 	write.table(paste(vertex_dir, "Ensembl_gene.bed", sep=""), sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)