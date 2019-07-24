library(tidyverse)

gwava_result_path <- './experiment_result/gwava_performance_CERENKOV2_1337.tsv'
c1_result_path <- './experiment_result/c1_cross_validate_xv_report.tsv'
c2_result_path <- './experiment_result/c2_cross_validate_xv_report.tsv'  # C1 + LS
c3_result_path <- './experiment_result/c3_cross_validate_xv_report.tsv'  # C1 + LS + N2V

gwava_result <- read_tsv(gwava_result_path)
c1_result <- read_tsv(c1_result_path)
c2_result <- read_tsv(c2_result_path)
c3_result <- read_tsv(c3_result_path)

# OMG "repeat" is a keyword in R (for Repeat Loops)
# Use df$`repeat` instead of df$repeat 
# I'll just rename this column
gwava_result <- gwava_result %>% rename(rep = `repeat`) 
c1_result <- c1_result %>% rename(rep = `repeat`) 
c2_result <- c2_result %>% rename(rep = `repeat`) 
c3_result <- c3_result %>% rename(rep = `repeat`)

gwava_result <- gwava_result %>% add_column(feat_set = "GWAVA") 
c1_result <- c1_result %>% add_column(feat_set = "CERENKOV") 
c2_result <- c2_result %>% add_column(feat_set = "CERENKOV+LS") 
c3_result <- c3_result %>% add_column(feat_set = "CERENKOV3") 

cv_result <- bind_rows(c1_result, c2_result, c3_result)

# AVGRANK is a the-lower-the-better metric
# 	and scikit-learn will automatically make it negative. 
# Flipping the sign will restore their values.
cv_result$test_AVGRANK <- -1 * cv_result$test_AVGRANK
cv_result$train_AVGRANK <- -1 * cv_result$train_AVGRANK

# GWAVA results are from R and the AVGRANKs are positive already
# No need to flip the sign 
cv_result <- bind_rows(gwava_result, cv_result)

# Ensure the order on x-axis
cv_result$feat_set <- factor(cv_result$feat_set, 
							 levels=c("GWAVA", "CERENKOV", "CERENKOV+LS", "CERENKOV3"))

font_size <- 9
out_dir <- "figures/"

##### AUPRC #####

output_file_name <- "AUPRC_errorbar"

p_auprc <- ggplot(cv_result, aes(x=feat_set, y=test_AUPRC)) + 
	# geom_boxplot() + 
	stat_summary(fun.y=mean, geom="point", size=2) + 
	stat_summary(fun.data = mean_se, fun.args = list(mult = 2), geom = "errorbar", width=0.5) + 
	ylab("AUPRC") + 
	theme_classic() + 
	theme(axis.title.x=element_blank(), 
		  axis.text.x=element_text(size=font_size, angle=45, hjust=1, vjust=1), 
		  axis.text.y=element_text(size=font_size))

ggsave(paste0(out_dir, output_file_name, ".pdf"), plot=p_auprc, width=2, height=3)
ggsave(paste0(out_dir, output_file_name, ".png"), plot=p_auprc, width=2, height=3)

##### AUROC #####

output_file_name <- "AUROC_errorbar"

p_auroc <- ggplot(cv_result, aes(x=feat_set, y=test_AUROC)) + 
	# geom_boxplot() + 
	stat_summary(fun.y=mean, geom="point", size=2) + 
	stat_summary(fun.data = mean_se, fun.args = list(mult = 2), geom = "errorbar", width=0.5) + 
	ylab("AUROC") + 
	theme_classic() + 
	theme(axis.title.x=element_blank(), 
		  axis.text.x=element_text(size=font_size, angle=45, hjust=1, vjust=1), 
		  axis.text.y=element_text(size=font_size))

ggsave(paste0(out_dir, output_file_name, ".pdf"), plot=p_auroc, width=2, height=3)
ggsave(paste0(out_dir, output_file_name, ".png"), plot=p_auroc, width=2, height=3)

##### AVGRANK #####

output_file_name <- "AVGRANK_errorbar"

p_avgrank <- ggplot(cv_result, aes(x=feat_set, y=test_AVGRANK)) + 
	# geom_boxplot() + 
	stat_summary(fun.y=mean, geom="point", size=2) + 
	stat_summary(fun.data = mean_se, fun.args = list(mult = 2), geom = "errorbar", width=0.5) + 
	ylab("AVGRANK") + 
	theme_classic() + 
	theme(axis.title.x=element_blank(), 
		  axis.text.x=element_text(size=font_size, angle=45, hjust=1, vjust=1), 
		  axis.text.y=element_text(size=font_size))

ggsave(paste0(out_dir, output_file_name, ".pdf"), plot=p_avgrank, width=2, height=3)
ggsave(paste0(out_dir, output_file_name, ".png"), plot=p_avgrank, width=2, height=3)
