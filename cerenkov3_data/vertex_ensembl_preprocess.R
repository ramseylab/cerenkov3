# This script aims to:
# 	1. rename the columns of `mart_export.txt`
# 	2. since Ensembl uses 1-based coordinate system while we use 0-based, 
# 		revise `TSS`, `gene_start` and `transcript_start` columns
#	3. add "chr" prefix to "chromosome" column, 
# 		e.g. "1" => "chr1"
# 	4. save the revised dataframe to a TSV file

library(dplyr)
library(readr)
library(magrittr)

ensembl_dir <- "./vertex/gene/Ensembl/"

ensembl_gene_df <- read_tsv(paste(ensembl_dir, "mart_export.txt", sep=""), col_types=cols(
	`Gene stable ID` = col_character(),
	Strand = col_integer(),
	`Chromosome/scaffold name` = col_character(),
	`Transcription start site (TSS)` = col_integer(),
	`Transcript start (bp)` = col_integer(),
	`Transcript end (bp)` = col_integer(),
	`Gene start (bp)` = col_integer(),
	`Gene end (bp)` = col_integer()
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
# Otherwise use `col_double()` when reading the TSV
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
# Extra to TSS: e.g. chr1:100-100 => chr1:99-100
ensembl_gene_df <- ensembl_gene_df %>% mutate(TSS_start = TSS - 1, 
											  TSS_end = TSS)

ensembl_gene_df %>% 
	select(-TSS) %>%  # Delete original `TSS` column				
	write_tsv(paste(ensembl_dir, "p1_Ensembl.tsv", sep=""))  # Save the final dataframe