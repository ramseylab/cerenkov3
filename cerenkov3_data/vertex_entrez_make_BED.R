# This script aims to:
    # 1. use `p2_Ensembl_promoter.bed`, `p2_Ensembl_transcript.bed` and `p2_Ensembl_TSS.bed` as inputs
    # 2. map Ensembl IDs to Entrez IDs using `p2_Ensembl_x_Entrez.tsv`
	# 3. save corresponding Entrez BED files
# Note: BED format is not strictly followed here.

library(dplyr)
library(readr)
library(magrittr)
library(tidyr)

read_bed <- function(path) {
	read_tsv(path, col_names=c("chrom", "chromStart", "chromEnd", "name", "id", "strand"), 
			 col_types=cols(chrom = col_character(),
							chromStart = col_integer(),
							chromEnd = col_integer(),
							name = col_character(),
							id = col_character(),
							strand = col_integer()))
}

ensembl_dir <- "./vertex/gene/Ensembl/"
ensembl_prm_bed <- read_bed(paste(ensembl_dir, "p2_Ensembl_promoter.bed", sep=""))
# ensembl_tx_bed <- read_bed(paste(ensembl_dir, "p2_Ensembl_transcript.bed", sep=""))
ensembl_tss_bed <- read_bed(paste(ensembl_dir, "p2_Ensembl_TSS.bed", sep=""))

entrez_dir <- "./vertex/gene/Entrez/"
e2e_map <- read_tsv(paste(entrez_dir, "p2_Ensembl_x_Entrez.tsv", sep=""), 
					col_types=cols(Ensembl_Gene_ID = col_character(),
								   Entrez_Gene_ID = col_character()))

make_entrez_bed <- function(ensembl_entrez_map, ensembl_bed) {
	ensembl_entrez_map %>% inner_join(ensembl_bed, by=c("Ensembl_Gene_ID" = "id"), copy=TRUE) %>%
			               select(chrom, chromStart, chromEnd, name, Entrez_Gene_ID, strand) %>%
						   distinct(.keep_all = TRUE) %>% 
						   arrange(chrom, chromStart)
}

entrez_prm_bed <- make_entrez_bed(e2e_map, ensembl_prm_bed)
# entrez_tx_bed <- make_entrez_bed(e2e_map, ensembl_tx_bed)
entrez_tss_bed <- make_entrez_bed(e2e_map, ensembl_tss_bed)

entrez_prm_bed %>% 
	write.table(paste(entrez_dir, "p3_Entrez_promoter.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)
# entrez_tx_bed %>% 
# 	write.table(paste(entrez_dir, "p3_Entrez_transcript.bed", sep=""), sep = "\t", 
# 				quote = FALSE, row.names = FALSE, col.names = FALSE)
entrez_tss_bed %>% 
	write.table(paste(entrez_dir, "p3_Entrez_TSS.bed", sep=""), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

