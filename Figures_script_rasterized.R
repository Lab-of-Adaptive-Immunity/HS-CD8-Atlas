# Author: Juraj Michalik
# Date: 15/11/2021

# This script generates figures used in manuscript written by Tsyklauri et al., 2022.
# The license applicable to this script (as well as others in this directory) can be consulted in GitHub repository (MIT, see respective file for details).
#
# Some of data-sets used in this atlas are downloaded from 10X site and as such are subject to their license;
# please see the list of used data-sets and consult their license (Creative Commons License).

# generates plots for Oksana's paper
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(dbscan)
library(scatterpie)
library(reshape2)
library(ggrepel)
library(scds)
library(SingleCellExperiment)
library(stringr)
library(svglite)
library(ggrastr)
library(scales)
library(ggnewscale)

# load data (MAIT)
donors.integrated <- readRDS('Donors_Datasets/Donors_Integrated_with_MAIT_final.rds')

# select top 5 invariants
frequency.table <- donors.integrated@meta.data %>%
  select(MAIT.motif) %>% group_by(MAIT.motif) %>%
  summarise(count = n()) %>%
  filter(MAIT.motif != 'Other motif') %>%
  arrange(-count) %>% slice_head(n = 5)

# extract and reformat needed data

data.expressions <- GetAssayData(donors.integrated)
data.expressions <- as.data.frame(t(data.expressions[rownames(data.expressions) %in% c('CD8A', 'CD8B', 'KLRB1'),]))
data.expressions[data.expressions < 0] <- 0 # min.cutoff = 0

donors.kept.data <- 
  cbind(as.data.frame(donors.integrated@reductions$umap@cell.embeddings), 
        donors.integrated@meta.data[c('seurat_clusters', 'optics', 'MAIT.motif')], 
        data.expressions) %>%
  mutate(MAIT.motif = ifelse(MAIT.motif %in% frequency.table$MAIT.motif,
                             MAIT.motif, 'Other motif'),
         pt.size = ifelse(MAIT.motif == 'Other motif', 0.1, 0.3))

# make desired plots
n_col <- length(unique(donors.kept.data$seurat_clusters))

p_clust <- 
  ggplot(donors.kept.data, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = as.factor(seurat_clusters)), size = 0.1)) +
  ggtitle('') +
  scale_color_manual(values = c(hue_pal()(n_col))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

n_col <- length(unique(donors.kept.data$optics))
colpal <- c(hue_pal()(n_col))
colpal[1] <- 'gray80'
colpal[4] <- 'firebrick2' 

p_clust_opt <- 
  ggplot(donors.kept.data, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = as.factor(optics)), size = 0.1)) +
  ggtitle('OPTICS') +
  scale_color_manual(values = colpal) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_klrb1 <- 
  ggplot(donors.kept.data, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = KLRB1), size = 0.1)) +
  ggtitle('Klrb1') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

n_col <- length(unique(donors.kept.data$MAIT.motif)) - 1

p_MAIT_invariants <-
  ggplot(donors.kept.data, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = MAIT.motif), size = donors.kept.data$pt.size)) +
  ggtitle('CAV.DSNYQLIW (5 most frequent)') +
  scale_color_manual(values = c(rainbow(n_col), '#CCCCCC66')) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  

# do same for MAIT-less data

donors.integrated.maitless <- readRDS('Donors_Datasets/Donors_Integrated_no_MAIT_final.rds')

# extract and reformat needed data

data.expressions <- GetAssayData(donors.integrated.maitless)
data.expressions <- as.data.frame(t(data.expressions[rownames(data.expressions) %in% c('IFITM3', 'IL7R', 'KLRD1',
                                                                                       'SELL', 'CD27', 'CD28', 'CCR7',
                                                                                       'FAS', 'IL2RB', 'GZMA', 'GZMB', 
                                                                                       'GZMK', 'GNLY', 'CD7'),]))
data.expressions[data.expressions < 0] <- 0 # min.cutoff = 0

donors.kept.data.nm <- 
  cbind(as.data.frame(donors.integrated.maitless@reductions$umap@cell.embeddings), 
        donors.integrated.maitless@meta.data[c('seurat_clusters', 'optics')], 
        data.expressions)

n_col <- length(unique(donors.kept.data.nm$seurat_clusters))

p_clust_nm <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = as.factor(seurat_clusters)), size = 0.1)) +
  ggtitle('') +
  scale_color_manual(values = c(hue_pal()(n_col))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

n_col <- length(unique(donors.kept.data.nm$optics))
colpal <- c(hue_pal()(n_col))
colpal[1] <- 'gray80'
colpal[4] <- 'firebrick2'            

p_clust_nm_opt <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = as.factor(optics)), size = 0.1)) +
  ggtitle('OPTICS') +
  scale_color_manual(values = colpal) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
p_ifitm3 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = IFITM3), size = 0.1)) +
  ggtitle('Ifitm3') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_il7r <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = IL7R), size = 0.1)) +
  ggtitle('Il7r') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_klrd1 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = KLRD1), size = 0.1)) +
  ggtitle('Klrd1') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_sell <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = SELL), size = 0.1)) +
  ggtitle('Cd62l') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd27 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CD27), size = 0.1)) +
  ggtitle('Cd27') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd28 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CD28), size = 0.1)) +
  ggtitle('Cd28') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_ccr7 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CCR7), size = 0.1)) +
  ggtitle('Ccr7') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd95 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = FAS), size = 0.1)) +
  ggtitle('Cd95') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd122 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = IL2RB), size = 0.1)) +
  ggtitle('Cd122') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_gzma <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = GZMA), size = 0.1)) +
  ggtitle('Gzma') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_gzmb <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = GZMB), size = 0.1)) +
  ggtitle('Gzmb') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_gzmk <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = GZMK), size = 0.1)) +
  ggtitle('Gzmk') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_gnly <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = GNLY), size = 0.1)) +
  ggtitle('Gnly') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd7 <- 
  ggplot(donors.kept.data.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CD7), size = 0.1)) +
  ggtitle('Cd7') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# antibody plotting
donors.subset <- subset(donors.integrated.maitless, orig.ident %in% c("Donor1_HS_CD8","Donor2_HS_CD8","Donor3_HS_CD8","Donor4_HS_CD8",
                                                                      "Healthy_Donor4_3pe_PBMC_10k_2S", "Healthy_Donor_5pe_PBMC_10k",
                                                                      "Healthy_Donor3_5pe_PBMC_10k", "Healthy_Donor4_5pe_PBMC_1k",
                                                                      "Healthy_Donor5_5pe_PBMC_10k"))

data.atb <- GetAssayData(donors.subset, assay = 'Antibodies')
data.atb <- as.data.frame(t(data.atb[rownames(data.atb) %in% c('CD45RA', 'CD45RO', 'CD8a'),]))
data.atb[data.atb < 0] <- 0 # min.cutoff = 0

donors.kept.atb.nm <- 
  cbind(as.data.frame(donors.subset@reductions$umap@cell.embeddings), 
        donors.subset@meta.data[c('seurat_clusters', 'optics')], 
        data.atb)

p_cd45ra <- 
  ggplot(donors.kept.atb.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CD45RA), size = 0.1)) +
  ggtitle('CD45RA') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd45ro <- 
  ggplot(donors.kept.atb.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CD45RO), size = 0.1)) +
  ggtitle('CD45RO') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p_cd8a <- 
  ggplot(donors.kept.atb.nm, aes(x = UMAP_1, y = UMAP_2)) +
  ggrastr::rasterise(geom_point(aes(col = CD8a), size = 0.1)) +
  ggtitle('Cd8a') +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.title=element_blank(),
        legend.key=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# prepare data for plotting clones
donor.gradient.data <- cbind(donors.integrated.maitless@meta.data[c('orig.ident', 'clonotype.repeated', 'cdr3.any', 'cdr3A.any', 'cdr3B.any')],
                             donors.integrated.maitless@reductions$umap@cell.embeddings) %>%
  arrange(-row_number()) %>%
  group_by(clonotype.repeated)

donor.gradient.data.repeated <- donor.gradient.data %>% filter(clonotype.repeated != 'Other') %>% 
  mutate(clone.occurence = n())
donor.gradient.data.other <- donor.gradient.data %>%
  mutate(other.cells = case_when(cdr3A.any == 'CDR3A' & cdr3B.any == 'CDR3B' ~ 'Complete',
                                 cdr3A.any == 'CDR3A' & cdr3B.any != 'CDR3B' ~ 'CDR3A Only',
                                 cdr3A.any != 'CDR3A' & cdr3B.any == 'CDR3B' ~ 'CDR3B Only',
                                 cdr3A.any != 'CDR3A' & cdr3B.any != 'CDR3B' ~ 'No VDJ'),
         other.cells = ifelse(grepl('_3pe_', orig.ident), 'No VDJ avail.', other.cells))
donor.gradient.data.other <- donor.gradient.data %>%
  mutate(global.info = ifelse(cdr3A.any == 'CDR3A' & cdr3B.any == 'CDR3B', 'CDR3', 'No VDJ'))
donor.gradient.hasdata <- donor.gradient.data %>%
  filter(clonotype.repeated == 'Other' & (cdr3A.any == 'CDR3A' | cdr3B.any == 'CDR3B'))


# clone plot
clone_plot <- 
  ggplot() +
  ggrastr::rasterise(geom_point(data = donor.gradient.data.other,
                                aes(x = UMAP_1, y = UMAP_2,
                                    col = other.cells), size = 0.1)) +
  scale_color_manual(name = 'Cells with unique \nor incomplete VDJ',
                     values = c('salmon', 'seagreen3', 'orchid2', 'gray40', 'black'),
                     guide = guide_legend(override.aes = list(size=3))) +
  new_scale_color() +
  ggrastr::rasterise(geom_point(data = donor.gradient.data.repeated,
                                aes(x = UMAP_1, y = UMAP_2,
                                    col = clone.occurence), size = 0.1)) +
  scale_color_gradient(name = 'Repeating clone \nfrequency',
                       low = "salmon", high = "firebrick4") +
  ggtitle('Human Atlas VDJ') + 
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.key = element_blank(),
        legend.title = element_text(face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

clone_plot_oth <- 
  ggplot() +
  ggrastr::rasterise(geom_point(data = donor.gradient.data.other,
                                aes(x = UMAP_1, y = UMAP_2,
                                    col = global.info), size = 0.1)) +
  scale_color_manual(name = 'Cell VDJ \ninformation',
                     values = c('goldenrod3', 'black'),
                     labels = c(expression('CDR3' ~ alpha ~ '+ CDR3' ~ beta), 'Partial/absent'),
                     guide = guide_legend(override.aes = list(size=3))) +
  ggtitle('Cells from Human Atlas with complete VDJs') + 
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.key = element_blank(),
        legend.title = element_text(face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

clone_plot_rep <- 
  ggplot() +
  ggrastr::rasterise(geom_point(data = donor.gradient.hasdata,
                                aes(x = UMAP_1, y = UMAP_2,
                                    col = cdr3.any), size = 0.1)) +
  scale_color_manual(name = '',
                     values = c('gray50'),
                     labels = c('Other cells with\npartial/complete\nVDJ information'),
                     guide = guide_legend(override.aes = list(size=3))) +
  new_scale_color() +
  ggrastr::rasterise(geom_point(data = donor.gradient.data.repeated,
                                aes(x = UMAP_1, y = UMAP_2,
                                    col = clone.occurence), size = 0.1)) +
  scale_color_gradient(name = 'Clone \noccurence',
                       low = "salmon", high = "firebrick4") +
  ggtitle('Human Atlas VDJ - repeating clones') + 
  theme(plot.title = element_text(hjust = 0.5,  face = 'bold.italic'),
        legend.key = element_blank(),
        legend.title = element_text(face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# save everything as svg
ggsave(file="General/Atlas_1_with_MAIT_Clustering.svg", plot=p_clust, width=10, height=10)
ggsave(file="General/Atlas_1_with_MAIT_Clustering_Optics.svg", plot=p_clust_opt, width=10, height=10)
ggsave(file="General/Atlas_1_with_MAIT_Klrb1.svg", plot=p_klrb1, width=10, height=10)
ggsave(file="General/Atlas_1_with_MAIT_invariants.svg", plot=p_MAIT_invariants, width=10, height=10)

ggsave(file="General/Atlas_1_no_MAIT_Clustering.svg", plot=p_clust_nm , width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Clustering_Optics.svg", plot=p_clust_nm_opt, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Ifitm3.svg", plot=p_ifitm3, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Il17r.svg", plot=p_il7r, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Klrd1.svg", plot=p_klrd1, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Sell.svg", plot=p_sell, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd27.svg", plot=p_cd27, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd28.svg", plot=p_cd28, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd95.svg", plot=p_cd95, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Ccr7.svg", plot=p_ccr7, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd122.svg", plot=p_cd122, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Gzma.svg", plot=p_gzma, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Gzmb.svg", plot=p_gzmb, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Gzmk.svg", plot=p_gzmk, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Gnly.svg", plot=p_gnly, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd7.svg", plot=p_cd7, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd45ra.svg", plot=p_cd45ra, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd45ro.svg", plot=p_cd45ro, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Cd8a.svg", plot=p_cd8a, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Clones.svg", plot=clone_plot, width=11, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Clones_rep.svg", plot=clone_plot_rep, width=10, height=10)
ggsave(file="General/Atlas_1_no_MAIT_Clones_other.svg", plot=clone_plot_oth, width=10, height=10)