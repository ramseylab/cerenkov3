# This script aims to:
# 	1. use `4D_Genome_df.tsv` as input
# 	2. save "InteractorAChr", "InteractorAStart", "InteractorAEnd" and "InteractionID"
# 		to a BED file, `InteractorA.bed`
# 	3. save "InteractorBChr", "InteractorBStart", "InteractorBEnd" and "InteractionID"
# 		to a BED file, `InteractorB.bed`

library(dplyr)
library(readr)
library(magrittr)

nd_genome_df <- read_tsv("./edge/snp-gene/Chromatin-Interaction/4DGenome/4DGenome_HomoSapiens_hg19_preprocessed.tsv")

nd_genome_df %>% 
	select(InteractorAChr, InteractorAStart, InteractorAEnd, InteractionID) %>%
	write.table("./edge/snp-gene/Chromatin-Interaction/4DGenome/4DGenome_InteractorA.bed", sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

nd_genome_df %>% 
	select(InteractorBChr, InteractorBStart, InteractorBEnd, InteractionID) %>%
	write.table("./edge/snp-gene/Chromatin-Interaction/4DGenome/4DGenome_InteractorB.bed", sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)
