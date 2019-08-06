library(tidyverse)

feat_path <- './cerenkov3_data/vertex/SNP/osu19_cerenkov_feat_mat_plus_group_size.tsv'
feat <- read_tsv(feat_path)

feat <- feat %>% select(label, group_size)
feat$label <- as.factor(feat$label)

out_dir <- "figures/"

save_plot <- function(out_dir, output_file_name, plot_obj, width, height) {
	ggsave(paste0(out_dir, output_file_name, ".pdf"), plot=plot_obj, width=width, height=height, device="png")
	ggsave(paste0(out_dir, output_file_name, ".png"), plot=plot_obj, width=width, height=height, device="png")
	# ggsave(paste0(out_dir, output_file_name, ".eps"), plot=plot_obj, width=width, height=height, device=cairo_ps)
	ggsave(paste0(out_dir, output_file_name, ".tif"), plot=plot_obj, width=width, height=height, device="tiff")
}

##### Density Plots #####

output_file_name <- "LS_density"

p_density <- ggplot(feat, aes(x=group_size, color=label, linetype=label)) + 
	stat_density(geom="line", position="identity", alpha=0.8, show.legend=TRUE) + 
	scale_linetype_manual(breaks=c("0", "1"), 
						  labels=c("cSNP", "rSNP"),
						  values=c("longdash", "solid")) +
	scale_color_manual(breaks=c("0", "1"), 
					   labels=c("cSNP", "rSNP"),
					   values=c("tomato", "blue")) + 
	xlab("Locus Size") + 
	scale_y_continuous(limits = c(0,0.0175), expand = c(0, 0)) + # remove the gap between baseline and x-axis
	theme_classic() + theme(text = element_text(size=9), 
							axis.text.y = element_blank(), 
							axis.ticks.y = element_blank(), 
							legend.position = c(0.8, 0.5))

save_plot(out_dir, output_file_name, p_density, 3, 2)

