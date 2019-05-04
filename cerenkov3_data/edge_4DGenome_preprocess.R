# This script aims to:
# 	1. rename column "Cell/Tissue" => "Cell_Tissue"
# 	2. show that there is only one unique value, "hg19", in "Organism" column so we can ignore it
# 	3. remove interactor regions with HUGE sizes (due to possible errors from 4DGenome)
# 	4. add a new column of "InteractionID"
# 	5. save the revised dataframe to another TSV file
# 	6. save "InteractorAChr", "InteractorAStart", "InteractorAEnd" and "InteractionID" to a BED file, `InteractorA.bed` (for BedTools)
# 	7. save "InteractorBChr", "InteractorBStart", "InteractorBEnd" and "InteractionID" to a BED file, `InteractorB.bed` (for BedTools)

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(magrittr)
library(purrr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(stringr)

res_dir <- "./resource/4DGenome/"

df_4DGenome <- read_tsv(paste0(res_dir, "4DGenome_HomoSapiens_hg19.txt"), col_types=cols(
	InteractorAStart = col_integer(),
	InteractorAEnd = col_integer(),
	InteractorBStart = col_integer(),
	InteractorBEnd = col_integer(),
	Agene = col_character(),
	Bgene = col_character(),
	Detection_Method = col_character(),
	Confidence_Score1 = col_double(),
	Confidence_Score2 = col_double(),
	Contact_Frequency = col_integer(),
	Pubmed_ID = col_character()
))

df_4DGenome <- df_4DGenome %>% rename(Cell_Tissue = `Cell/Tissue`)

# There should be only one unique value in "Organism" column
stopifnot(df_4DGenome %>% distinct(Organism) %>% count == 1)
stopifnot(df_4DGenome %>% distinct(Organism) == "hg19")

# There are 4 rows with ridiculous InteractorA sizes
# 	Detection_Method	Pubmed_ID	InteractorAChr	InteractorAStart	InteractorAEnd	InteractorBChr	InteractorBStart	InteractorBEnd	InteractorALength
# 	Hi-C	24141950	chr2	9231530	92321964	chr2	92321965	92324000	83090434
# 	Hi-C	24141950	chr2	13308860	133025619	chr2	133025620	133030966	119716759
# 	Hi-C	24141950	chr2	13309013	133038761	chr2	133012117	133019012	119729748
# 	Hi-C	24141950	chr2	13302117	133038761	chr2	133008860	133012116	119736644
# All other interactors are shorter than 500000bp
interactor_size_thold <- 500000
df_4DGenome <- df_4DGenome %>% 
	mutate(InteractorASize = InteractorAEnd - InteractorAStart, 
		   InteractorBSize = InteractorBEnd - InteractorBStart) %>% 
	filter(InteractorASize <= interactor_size_thold & 
		   InteractorBSize <= interactor_size_thold)


# df_4DGenome %>% filter(!is.na(Agene)) %>% filter(!is.na(Bgene)) %>% select(Agene, Bgene) %T>% write_tsv("foo.tsv")

# Find trans-chrom interactors
# df_4DGenome %>% filter(InteractorAChr != InteractorBChr)

df_4DGenome <- df_4DGenome %>% 
	group_by(Pubmed_ID) %>% 
	mutate(inner_id = row_number()) %>%  # "inner_id" is the ID within a group
	unite(InteractionID, Pubmed_ID, inner_id, remove=TRUE) %>%  # remove "Pubmed_ID" and "inner_id" column
	ungroup %>%
	select(-Organism, -InteractorASize, -InteractorBSize) 
	
df_4DGenome	%>% 
	write_tsv(paste0(res_dir, "4DGenome_HomoSapiens_hg19_preprocessed.tsv"))

# Save BED files from the preprocessed data
df_4DGenome %>% 
	select(InteractorAChr, InteractorAStart, InteractorAEnd, InteractionID) %>%
	write.table(paste0(res_dir, "InteractorA.bed"), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)

df_4DGenome %>% 
	select(InteractorBChr, InteractorBStart, InteractorBEnd, InteractionID) %>%
	write.table(paste0(res_dir, "InteractorB.bed"), sep = "\t", 
				quote = FALSE, row.names = FALSE, col.names = FALSE)