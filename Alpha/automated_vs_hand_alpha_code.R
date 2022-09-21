rm(list=ls());if(is.null(dev.list()["RStudioGD"])){} else {dev.off(dev.list()["RStudioGD"])};cat("\014")
#_________________________________________________________________________________________________________________________
#Alpha Diversity - (3G-ESP/Hand) Method Comparison 
#_________________________________________________________________________________________________________________________
#Load dependencies 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomformat")
BiocManager::install("phyloseq")

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(tidyverse)
library(biomformat)
library(phyloseq)
library(qiime2R)
library(tibble)
library(dplyr)
library(ggpubr)

#Set working directory to the location of this R code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Create phyloseq file 
#Import ASV table (rarefied)
ASVtable.rare <-read_qza("rarefied_table2019_16s_Feb2022_filtered.qza")
ASVtable.rare <- ASVtable.rare$data

#Import taxonomy file, specify taxa headers, and remove confidence values
ASVtaxa <- read_qza("taxonomy.qza")

taxtable <- ASVtaxa$data %>% separate(Taxon, sep=";",
                                      c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
taxtable$Confidence <- NULL

#Read in metadata file and convert empty cells to "NA"
sample_info_table <- read.table("metadata_2019.tsv",  sep ='\t', header = TRUE,
                                row.names = 1, na.strings = c("", "NA"))

#remove everything except Lake Erie samples
sample_info_table <- subset(sample_info_table, description == "Lake Erie sample" | description == "Lake Erie hand sample")

physeq.rare <- phyloseq(otu_table(ASVtable.rare, taxa_are_rows= T), 
                        tax_table(as.data.frame(taxtable) %>% column_to_rownames("Feature.ID") %>% 
                                    as.matrix()), sample_data(sample_info_table))

#Calculate shannon diversity for 16S sequences
ESP_hand_16S_shannon <- estimate_richness(physeq.rare, split = TRUE, measures = c("Shannon"))
colnames(ESP_hand_16S_shannon)[colnames(ESP_hand_16S_shannon) == "Shannon"] <- "shannon_16s"
ESP_hand_16S_shannon %>% tibble::rownames_to_column(.,"sample_name") -> ESP_hand_16S_shannon

#Load calculated nonpareil diversity (Nd) for metagenomic sequences
nonpareil_2019_Nd <- read.csv("nonpareil_output.csv", header = TRUE,
                              row.names = 1, na.strings = c("", "NA"))

##Edit nonpareil sample names to find a common name with metadata
row.names(nonpareil_2019_Nd) <- sub("_R1_trim", "", row.names(nonpareil_2019_Nd))
colnames(nonpareil_2019_Nd)[colnames(nonpareil_2019_Nd) == "diversity"] <- "Nd"
nonpareil_2019_Nd <- tibble::rownames_to_column(nonpareil_2019_Nd, "sample_id_shotgun_short")

##Extract metadata, mates only
WLE_3GESP19_mates <- sample_info_table %>% filter(!is.na(mate_label))
WLE_3GESP19_mates <- rownames_to_column(WLE_3GESP19_mates, "sample_name")
WLE_3GESP19_mates <- WLE_3GESP19_mates %>% add_column(sample_id_shotgun_short = gsub("_........-........", "", WLE_3GESP19_mates$sample_id_shotgun))
#Add/adjust metadata
ESP_hand_alpha_comb <- right_join(ESP_hand_16S_shannon,WLE_3GESP19_mates,by = "sample_name")
ESP_hand_alpha_Nd_comb <- right_join(nonpareil_2019_Nd,ESP_hand_alpha_comb,by = "sample_id_shotgun_short")

ESP_hand_alpha_Nd_comb$filter_size = str_sub(ESP_hand_alpha_Nd_comb$sample_name,-2) 
ESP_hand_alpha_Nd_comb_sub <- ESP_hand_alpha_Nd_comb[,c("sample_name","shannon_16s", "Nd", "description", "sample_type", "filter_size")]

#_________________________________________________________________________________________________________________________
#Alpha Diversity Comparison - 16S amplicon 5.0um fraction
#_________________________________________________________________________________________________________________________
#5.0um samples
wilcox_test_16S_5.0 <- wilcox.test(subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Archive 5.0um")$shannon_16s, 
                                   subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Hand-sample 5.0um")$shannon_16s)
wilcox_test_16S_5.0 

#Wilcoxon rank sum exact test
#data:  subset(ESP_hand_16S_shannon_sub, sample_type == "Archive 5.0um")$Shannon and subset(ESP_hand_16S_shannon_sub, sample_type == "Hand-sample 5.0um")$Shannon
#W = 55, p-value = 0.7394

#_________________________________________________________________________________________________________________________
#Alpha Diversity Comparison - 16S amplicon 0.22um fraction
#_________________________________________________________________________________________________________________________
#0.22um samples
wilcox_test_16S_0.22 <- wilcox.test(subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Archive 0.22um")$shannon_16s, 
                                    subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Hand-sample 0.22um")$shannon_16s)
wilcox_test_16S_0.22 

#Wilcoxon rank sum exact test
#data:  subset(ESP_hand_16S_shannon_sub, sample_type == "Archive 0.22um")$Shannon and subset(ESP_hand_16S_shannon_sub, sample_type == "Hand-sample 0.22um")$Shannon
#W = 46, p-value = 0.7959

#_________________________________________________________________________________________________________________________
#Alpha Diversity Comparison - Nonpareil - Metagenomic 5.0um fraction
#_________________________________________________________________________________________________________________________
#5.0um samples 
wilcox_test_shotgun_50 <- wilcox.test(subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Archive 5.0um")$Nd, 
                                      subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Hand-sample 5.0um")$Nd)
wilcox_test_shotgun_50 
#Wilcoxon rank sum exact test

#data:  subset(ESP_hand_Nd_comb_50, sample_type == "Archive 5.0um")$diversity and subset(ESP_hand_Nd_comb_50, sample_type == "Hand-sample 5.0um")$diversity
#W = 34, p-value = 0.2475

#_________________________________________________________________________________________________________________________
#Alpha Diversity Comparison - Nonpareil - Metagenomic 0.22um fraction
#_________________________________________________________________________________________________________________________
#5.0um samples
wilcox_test_shotgun_22 <- wilcox.test(subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Archive 0.22um")$Nd, 
                                      subset(ESP_hand_alpha_Nd_comb_sub, sample_type == "Hand-sample 0.22um")$Nd)
wilcox_test_shotgun_22
#Wilcoxon rank sum exact test

#data:  subset(ESP_hand__shannon_comb_22, sample_type == "Archive 0.22um")$diversity and subset(ESP_hand__shannon_comb_22, sample_type == "Hand-sample 0.22um")$diversity
#W = 68, p-value = 0.1903

#_________________________________________________________________________________________________________________________
#Violin Plots---16S amplicon 5.0um and 0.22um fractions
#_________________________________________________________________________________________________________________________
#Change Lake Erie hand sample to "Hand sample" and Lake Erie sample to "3G-ESP"
ESP_hand_alpha_Nd_comb %>%
  mutate(description = ifelse(description == "Lake Erie hand sample", "Hand sample", "3G-ESP")) -> ESP_hand_alpha_Nd_comb

#Designate colors to sampling methods
Hand_ESP_colors = c("Hand sample" = "#6a3d9a", "3G-ESP" = "#33a02c")

#5.0um - 16s
data_16s_50 = subset(ESP_hand_alpha_Nd_comb, filter_size == 50)
plot_16s_50 <- ggplot(data_16s_50 , aes(x = description, y = shannon_16s, colour=description)) + 
  geom_point(size=1, position=position_dodge(width=1)) +
  geom_violin(alpha=0, position=position_dodge(width=1)) + 
  scale_color_manual(values = Hand_ESP_colors) + 
  ylim(0,7) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.line.y = element_line(size = ifelse(50 %in% c(2,4), 0, 0.5), color = "black"),
        axis.ticks.length.x = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.ticks.length.y = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.text.y = element_text(size = ifelse(50 %in% c(2,4), 0, 15)),
        axis.text.x = element_text(size = ifelse(50 %in% c(1,2), 0, 15)),
        axis.title = element_blank()) +
  annotate(geom = "text", label = paste0("16S Amplicon Sequencing\n5.0µm Filter", ": ", "p > 0.05"), x = 2.4, y = 6.8, size = 5, hjust = 1)+
  labs(y ="Shannon Index")

#0.22um - 16s
data_16s_22 = subset(ESP_hand_alpha_Nd_comb, filter_size == 22)
plot_16s_22 <- ggplot(data_16s_22, aes(x = description, y = shannon_16s, colour=description)) + 
  geom_point(size=1, position=position_dodge(width=1)) +
  geom_violin(alpha=0, position=position_dodge(width=1)) + 
  scale_color_manual(values = Hand_ESP_colors) + 
  ylim(0,7) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.line.y = element_line(size = ifelse(50 %in% c(2,4), 0, 0.5), color = "black"),
        axis.ticks.length.x = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.ticks.length.y = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.text.y = element_text(size = ifelse(50 %in% c(2,4), 0, 15)),
        axis.text.x = element_text(size = ifelse(50 %in% c(1,2), 0, 15)),
        axis.title = element_blank()) +
  annotate(geom = "text", label = paste0("16S Amplicon Sequencing\n0.22µm Filter", ": ", "p > 0.05"), x = 2.4, y = 6.8, size = 5, hjust = 1)+
  labs(y ="Shannon Index")

#_________________________________________________________________________________________________________________________
#Violin Plots---Metagenomic (shotgun) 5.0um and 0.22um fractions
#_________________________________________________________________________________________________________________________
#5.0um - shotgun
data_shotgun_50 = subset(ESP_hand_alpha_Nd_comb, filter_size == 50)
plot_shotgun_Nd_50 <- ggplot(data_shotgun_50, aes(x = description, y = Nd, colour=description)) + 
  geom_point(size=1, position=position_dodge(width=1)) +
  geom_violin(alpha=0, position=position_dodge(width=1)) + 
  scale_color_manual(values = Hand_ESP_colors) + 
  ylim(16,24) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.line.y = element_line(size = ifelse(50 %in% c(2,4), 0, 0.5), color = "black"),
        axis.ticks.length.x = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.ticks.length.y = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.text.y = element_text(size = ifelse(50 %in% c(2,4), 0, 15)),
        axis.text.x = element_text(size = ifelse(50 %in% c(1,2), 0, 15)),
        axis.title = element_blank()) +
  annotate(geom = "text", label = paste0("Metagenomic Sequencing\n5.0µm Filter", ": ", "p > 0.05"), x = 2.4, y = 23, size = 5, hjust = 1)+
  labs(y ="Nonpareil diversity (Nd)")

#0.22um - shotgun
data_shotgun_22 = subset(ESP_hand_alpha_Nd_comb, filter_size == 22)
plot_shotgun_Nd_22 <- ggplot(data_shotgun_22, aes(x = description, y = Nd, colour=description)) + 
  geom_point(size=1, position=position_dodge(width=1)) +
  geom_violin(alpha=0, position=position_dodge(width=1)) + 
  scale_color_manual(values = Hand_ESP_colors) + 
  ylim(16,24) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.line.y = element_line(size = ifelse(50 %in% c(2,4), 0, 0.5), color = "black"),
        axis.ticks.length.x = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.ticks.length.y = unit(ifelse(50 %in% c(2,4), 0, 0.2),"cm"),
        axis.text.y = element_text(size = ifelse(50 %in% c(2,4), 0, 15)),
        axis.text.x = element_text(size = ifelse(50 %in% c(1,2), 0, 15)),
        axis.title = element_blank()) +
  annotate(geom = "text", label = paste0("Metagenomic Sequencing\n0.22µm Filter", ": ", "p > 0.05"), x = 2.4, y = 23, size = 5, hjust = 1)+
  labs(y ="Nonpareil diversity (Nd)")

#Arrange all plots
#16s amplicon
plot_16s_amplicon_combined <- annotate_figure(ggpubr::ggarrange(plot_16s_50, plot_16s_22,
                                                           ncol = 2, nrow = 1,labels = c("(a)", "(b)"), 
                                                           label.x = 0.1, label.y = 0.985,
                                                           font.label = list(size = 20, color = "black", face = "plain")),
                                         # annotate_figure
                                         left = text_grob("Shannon Index", color = "black", size = 15, rot = 90),
                                         bottom = text_grob("Sampling Method", size = 15))
ggsave(plot_16s_amplicon_combined , height = 4, width = 8, path = ".", filename = "plot_16s_amplicon_combined.pdf", device = "pdf")

#Shotgun
plot_shotgun_combined <- annotate_figure(ggpubr::ggarrange(plot_shotgun_Nd_50, plot_shotgun_Nd_22,
                                                      ncol = 2, nrow = 1,labels = c("(c)", "(d)"), 
                                                      label.x = 0.1, label.y = 0.985,
                                                      font.label = list(size = 20, color = "black", face = "plain")),
                                    # annotate_figure
                                    left = text_grob("Nonpareil diversity (Nd)", color = "black", size = 15, rot = 90),
                                    bottom = text_grob("Sampling Method", size = 15))
ggsave(plot_shotgun_combined, height = 4, width = 8, path = ".", filename = "plot_shotgun_combined.pdf", device = "pdf")




