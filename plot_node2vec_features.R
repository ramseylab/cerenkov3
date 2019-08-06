library(tidyverse)

n2v_feat_path <- './cerenkov3_classifier/N2V_FEAT/Sgn_3.0_0.3_0.3_0.1_0.1_0.3_0.3_3.0_sum_2_6_12_4_4_8_True.tsv'
n2v_feat <- read_tsv(n2v_feat_path)

n2v_feat$label <- as.factor(n2v_feat$label)

out_dir <- "figures/"

save_plot <- function(out_dir, output_file_name, plot_obj, width, height) {
	ggsave(paste0(out_dir, output_file_name, ".pdf"), plot=plot_obj, width=width, height=height, device="png")
	ggsave(paste0(out_dir, output_file_name, ".png"), plot=plot_obj, width=width, height=height, device="png")
	# ggsave(paste0(out_dir, output_file_name, ".eps"), plot=plot_obj, width=width, height=height, device=cairo_ps)
	ggsave(paste0(out_dir, output_file_name, ".tif"), plot=plot_obj, width=width, height=height, device="tiff")
}

##### Scatter Plot #####

output_file_name <- "N2V_feat_scatter"

p_scatter <- ggplot(n2v_feat, aes(x=node2vec_1, y=node2vec_2, color=label)) + 
	geom_point(shape=18, alpha=0.2) + 
	scale_color_manual(breaks=c("0", "1"), 
					   labels=c("cSNP", "rSNP"),
					   values=c("tomato", "blue")) + 
	xlab("1st Node2vec Feature") + ylab("2nd Node2vec Feature") + 
	theme_classic() + theme(text = element_text(size=9))

save_plot(out_dir, output_file_name, p_scatter, 3, 2)

##### Density Plots #####

output_file_name <- "N2V_feat1_density"

p_density_1 <- ggplot(n2v_feat, aes(x=node2vec_1, color=label, linetype=label)) + 
	stat_density(geom="line", position="identity", alpha=0.8, show.legend=TRUE) + 
	scale_linetype_manual(breaks=c("0", "1"), 
	                      labels=c("cSNP", "rSNP"),
					      values=c("longdash", "solid")) + 
	scale_color_manual(breaks=c("0", "1"), 
	                   labels=c("cSNP", "rSNP"),
					   values=c("tomato", "blue")) + 
	xlab("1st Node2vec Feature") + 
	scale_y_continuous(limits = c(0,1.75), expand = c(0, 0)) + # remove the gap between baseline and x-axis
	theme_classic() + theme(text = element_text(size=9), 
							axis.text.y = element_blank(), 
							axis.ticks.y = element_blank(), 
							legend.position = c(0.8, 0.5))

save_plot(out_dir, output_file_name, p_density_1, 3, 2)


output_file_name <- "N2V_feat2_density"

p_density_2 <- ggplot(n2v_feat, aes(x=node2vec_2, color=label, linetype=label)) + 
	stat_density(geom="line", position="identity", alpha=0.8, show.legend=TRUE) + 
	scale_linetype_manual(breaks=c("0", "1"), 
	                      labels=c("cSNP", "rSNP"),
					      values=c("longdash", "solid")) + 
	scale_color_manual(breaks=c("0", "1"), 
	                   labels=c("cSNP", "rSNP"),
					   values=c("tomato", "blue")) + 
	xlab("2nd Node2vec Feature") + 
	scale_y_continuous(limits = c(0,1.75), expand = c(0, 0)) + # remove the gap between baseline and x-axis
	theme_classic() + theme(text = element_text(size=9), 
							axis.text.y = element_blank(), 
							axis.ticks.y = element_blank(), 
							legend.position = c(0.8, 0.5))

save_plot(out_dir, output_file_name, p_density_2, 3, 2)

## Put 2 density plots together

melt_n2v_feat <- n2v_feat %>% select(name, label, 
									 `1st Node2vec Feature`=node2vec_1, 
									 `2nd Node2vec Feature`=node2vec_2) %>% 
							  gather(`1st Node2vec Feature`, `2nd Node2vec Feature`, key = "feat_name", value = "feat_value")

output_file_name <- "N2V_feat_densities"

p_densities <- ggplot(melt_n2v_feat, aes(x=feat_value, color=label, linetype=label)) + 
	facet_wrap(~feat_name, nrow = 1, scales = "free_x", strip.position = "bottom") + 
	stat_density(geom="line", position="identity", alpha=0.8, show.legend=TRUE) + 
	scale_linetype_manual(breaks=c("0", "1"), 
	                      labels=c("cSNP", "rSNP"),
					      values=c("longdash", "solid")) + 
	scale_color_manual(breaks=c("0", "1"), 
	                   labels=c("cSNP", "rSNP"),
					   values=c("tomato", "blue")) + 
	scale_y_continuous(limits = c(0, 1.75), expand = c(0, 0)) + # remove the gap between baseline and x-axis
	theme_classic() + theme(text = element_text(size=9), 
							axis.title.x = element_blank(), 
							axis.text.y = element_blank(), 
							axis.ticks.y = element_blank(), 
							strip.text = element_text(size=9, vjust=2.5),
							strip.background = element_blank(), 
							strip.placement = "outside", 
							legend.position = c(0.4, 0.5))

save_plot(out_dir, output_file_name, p_densities, 6, 2)

