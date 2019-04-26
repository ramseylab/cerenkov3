# This script aims to:
# 	1. rename column "Cell/Tissue" => "Cell_Tissue"
# 	2. show that there is only one unique value, "hg19", in "Organism" column
# 		so we can ignore it
# 	3. remove interactor regions with HUGE sizes (may be errors from 4DGenome)
# 	4. add a new column of "InteractionID"
# 	5. save the revised dataframe to another TSV file

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(magrittr)
library(purrr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(stringr)

# well, a variable name cannot starts with a number
nd_genome_df <- read_tsv("./edge/snp-gene/Chromatin-Interaction/4DGenome/4DGenome_HomoSapiens_hg19.txt", col_types=cols(
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

nd_genome_df <- nd_genome_df %>% rename(Cell_Tissue = `Cell/Tissue`)

# There should be only one unique value in "Organism" column
stopifnot(nd_genome_df %>% distinct(Organism) %>% count == 1)
stopifnot(nd_genome_df %>% distinct(Organism) == "hg19")

# There are 4 rows with ridiculous InteractorA sizes
# 	Detection_Method	Pubmed_ID	InteractorAChr	InteractorAStart	InteractorAEnd	InteractorBChr	InteractorBStart	InteractorBEnd	InteractorALength
# 	Hi-C	24141950	chr2	9231530	92321964	chr2	92321965	92324000	83090434
# 	Hi-C	24141950	chr2	13308860	133025619	chr2	133025620	133030966	119716759
# 	Hi-C	24141950	chr2	13309013	133038761	chr2	133012117	133019012	119729748
# 	Hi-C	24141950	chr2	13302117	133038761	chr2	133008860	133012116	119736644
# All other interactors are shorter than 500000bp
interactor_size_thold <- 500000
nd_genome_df <- nd_genome_df %>% 
	mutate(InteractorASize = InteractorAEnd - InteractorAStart, 
		   InteractorBSize = InteractorBEnd - InteractorBStart) %>% 
	filter(InteractorASize <= interactor_size_thold & 
		   InteractorBSize <= interactor_size_thold)
									

# Does not work! More have more than 2 ","-separated values
# nd_genome_df %>% separate(Agene, sep=",", into=c("Agene_Name", "Agene_Id"))

# nd_genome_df %>% filter(!is.na(Agene)) %>% filter(!is.na(Bgene)) %>% select(Agene, Bgene) %T>% write_tsv("foo.tsv")

# Find trans-chrom interactors
# nd_genome_df %>% filter(InteractorAChr != InteractorBChr)

nd_genome_df %>% 
	group_by(Pubmed_ID) %>% 
	mutate(inner_id = row_number()) %>%  # "inner_id" is the ID within a group
	unite(InteractionID, Pubmed_ID, inner_id, remove=TRUE) %>%  # remove "Pubmed_ID" and "inner_id" column
	ungroup %>%
	select(-Organism, -InteractorASize, -InteractorBSize) %>% 
	write_tsv("./edge/snp-gene/Chromatin-Interaction/4DGenome/4DGenome_HomoSapiens_hg19_preprocessed.tsv")
