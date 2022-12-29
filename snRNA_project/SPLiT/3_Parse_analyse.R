library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pathfindR)
##################################################
# Perform DGE analysis between different a-syn bins, run pathfindR on them and generate the relevant plots


pbmc.all.da <- readRDS("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/pbmc_parse_da_subset.rds")
Idents(pbmc.all.da) <- pbmc.all.da$ASYN_range
table(Idents(pbmc.all.da))
table(pbmc.all.da@active.ident)
Idents(pbmc.all.da) <- pbmc.all.da$ASYN_range
table(Idents(pbmc.all.da))
pbmc.all.da@active.ident <- factor(pbmc.all.da@active.ident, 
                                   levels=c("0", "1", "2-3", "4-8", "9-15", "16-93"))


counts_s4 <- pbmc.all.da@assays$originalexp@counts
counts <- as.data.frame(counts_s4)
counts[1:5,1:5]
write.csv(counts, "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/da_subtype_counts.csv")

dge_p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/parse_var_genes_da_further_subset/"

# 0 1 2-3 9-93  4-8 
library(ggplot2)

pbmc.all.da.markers_all <- FindAllMarkers(pbmc.all.da)
head(pbmc.all.da.markers_all, 40) -> top40
# pbmc.all.da.markers_1_2_3 %>% top_n(-40, p_val_adj) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("Find all markers function")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("Find all markers function")



pbmc.all.da.markers_0_1 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.da.markers_0_1, 10)
write.csv(pbmc.all.da.markers_0_1, paste0(dge_p, "markers_0_vs_1.csv", ))

head(pbmc.all.da.markers_0_1, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 0 vs 1")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 0 vs 1")

###############################################################
asyn_0_vs_1_dotplot
pbmc.all.da.markers_1_2_3 <- FindMarkers(pbmc.all.da, ident.1 = "1", ident.2 = "2-3", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
write.csv(pbmc.all.da.markers_1_2_3, paste0(dge_p, "markers_1_vs_2_3.csv"))


head(pbmc.all.da.markers_1_2_3, 40) -> top40
# pbmc.all.da.markers_1_2_3 %>% top_n(-40, p_val_adj) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 1 vs 2-3")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 1 vs 2-3")

###############################################################
pbmc.all.da.markers_2_3__4_8 <- FindMarkers(pbmc.all.da, ident.1 = "2-3", ident.2 = "4-8", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
write.csv(pbmc.all.da.markers_2_3__4_8, paste0(dge_p, "markers_2_3_vs_4_8.csv"))

head(pbmc.all.da.markers_2_3__4_8, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 2-3 vs 4-8")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 2-3 vs 4-8")

###############################################################
table(Idents(pbmc.all.da))
pbmc.all.da.markers_4_8__16_93 <- FindMarkers(pbmc.all.da, ident.1 = "4-8", ident.2 = "16-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
write.csv(pbmc.all.da.markers_4_8__16_93, paste0(dge_p, "markers_4_8_vs_9_93.csv"))

head(pbmc.all.da.markers_4_8__9_93, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()

###############################################################
table(Idents(pbmc.all.da))
table(pbmc.all.da$ASYN_range)
Idents(pbmc.all.da) <- pbmc.all.da$ASYN_range

pbmc.all.da.markers_0__9_93 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "16-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
pbmc.all.da.markers_9_93_0 <- FindMarkers(pbmc.all.da, ident.1 = "16-93", ident.2 = "0", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)

head(pbmc.all.da.markers_0__9_93)
head(pbmc.all.da.markers_9_93_0)
#############################################################
pbmc.all.da.markers_0_1
table(Idents(pbmc.all.da))
pbmc.all.da.markers_0_1 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
pbmc.all.da.markers_0_2_3 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "2-3", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
pbmc.all.da.markers_0_4_8 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "4-8", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
pbmc.all.da.markers_0_9_15 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "9-15", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)

pbmc.all.da.markers_9_15__16_93 <- FindMarkers(pbmc.all.da, ident.1 = "9-15", ident.2 = "16-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
pbmc.all.da.markers_all <- FindAllMarkers(pbmc.all.da)
#############################################################

logFC <- as.vector(pbmc.all.da.markers_0__9_93$avg_log2FC)
FDR_p <- as.vector(pbmc.all.da.markers_0__9_93$p_val)
Gene_symbol <- as.vector(rownames(pbmc.all.da.markers_0__9_93))

pathfindr_df <- function(marker_df){
  logFC <- as.vector(marker_df$avg_log2FC)
  FDR_p <- as.vector(marker_df$p_val)
  Gene_symbol <- as.vector(rownames(marker_df))
  input_markers <- data.frame(Gene_symbol, logFC, FDR_p)
  return(input_markers)
}

input_markers_da_0_vs_93 <- pathfindr_df(pbmc.all.da.markers_0__9_93)

input_markers_da_0_vs_1 <- pathfindr_df(pbmc.all.da.markers_0_1)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_0_2_3)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_0_4_8)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_0__9_93)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_0_9_15)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_9_15__16_93)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_all)
input_markers_da <- pathfindr_df(pbmc.all.da.markers_9_93_0)

Enrichment_dotplot_da_asyn_0_vs_4_8

head(input_markers_da_0_vs_93)
library(pathfindR)
save_path = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/pathfindR_results"
output_markers_da_0_vs_1 <- run_pathfindR(input_markers_da_0_vs_1, output_dir = save_path, gene_sets = "KEGG")

output_markers <- run_pathfindR(input_markers_da, gene_sets = "KEGG")

# output_markers_da_0_vs_93 <- run_pathfindR(input_markers_da_0_vs_93, output_dir = save_path, gene_sets = "KEGG")


output_markers_da_0_vs_93
term_gene_heatmap(output_markers)
term_gene_graph(output_markers)
# UpSet_plot(output_markers)

# default settings
clustered_df <- cluster_enriched_terms(output_markers_da_0_vs_93)

# display the heatmap of hierarchical clustering
clustered_df <- cluster_enriched_terms(output_markers_da_0_vs_93, plot_hmap = TRUE)

# display the dendrogram and automatically-determined clusters
clustered_df <- cluster_enriched_terms(output_markers_da_0_vs_93, plot_dend = TRUE)

# change agglomeration method (default = "average") for hierarchical clustering
clustered_df <- cluster_enriched_terms(output_markers_da_0_vs_93, clu_method = "centroid")
###################################################
output_markers_da_0_vs_93[1,1]
plot_aav_load


write.csv(pbmc.all.da.markers_0__9_93, paste0(dge_p, "markers_0_vs_9_93.csv"))


head(pbmc.all.da.markers_0__9_93, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 0 vs 16-93")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 0 vs 16-93")


###########################################################
pbmc.all.da.markers_1__9_93 <- FindMarkers(pbmc.all.da, ident.1 = "1", ident.2 = "16-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
write.csv(pbmc.all.da.markers_1__9_93, paste0(dge_p, "markers_1_vs_9_93.csv"))

head(pbmc.all.da.markers_1__9_93, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 1 vs 16-93")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 1 vs 16-93")
###########################################################



pbmc.all.da.markers_0_41.93 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "41-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.da.markers_0_41.93, 40) <- top40
write.csv(pbmc.all.da.markers_0_41.93, paste0(dge_p, "markers_0_vs_41-93.csv"))
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()

#############################################################
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("ReactomePA")
# Gene expression analysis plots
NormalizeData(pbmc.all.da)
norm_data(data, norm.fun = "CNF", parameter.list = NULL, data.trans = "none")

counts_s4 <- pbmc.all.da@assays$originalexp@counts
counts <- as.data.frame(counts_s4)
# counts <- log2(counts)

counts < log(counts+1, base = 2)


counts_t = t(counts)
counts_t$ASYN_BINS <- pbmc.all.da$ASYN_range
counts_t$ASYN_BINS
colnames(counts_t)
counts_t["AAV-ASYN"]

plot_aav_load <- function(counts_t, aav, gene){
  gene_list = as.vector(unlist(colnames(counts_t)))
  
  Asyn_load <- as.data.frame(counts_t[, c(which(gene_list==c(aav)))])
  colnames(Asyn_load) <- aav
  Asyn_load_log <- log(Asyn_load+1, base = 2)
  
  Gene_df <- as.data.frame(counts_t[, c(which(gene_list==c(gene)))])
  colnames(Gene_df) <- gene
  Gene_df_log <- log(Gene_df+1, base = 2)
  
  asyn_x_gene <- cbind(Asyn_load_log, Gene_df_log)
  
  plot(asyn_x_gene[,1], asyn_x_gene[,2], main="gene expression relative to asyn load", xlab="Asyn load", ylab="Gene expression", pch=19) 
  
  abline(lm(asyn_x_gene[,2] ~ asyn_x_gene[,1]), col="red") # regression line (y~x) 
  
}

table(pbmc.all.da$Annotations_specific)
plot_aav_load(counts_t, "AAV-ASYN", "SOX6")
Idents(pbmc.all.da) <- pbmc.all.da$Annotations_specific

DotPlot(pbmc.all.da, features=c('TH', 'SLC6A3', 'CALB1', 'CALB2'))
Idents(pbmc.all.da) <- pbmc.all.da$ASYN_range
# NG2 not found
DotPlot(pbmc.all.da, features=c('RGS5', 'SIRT5', 'SIRT1', 'SIRT3', 'SPHK1', 'SOD2', 'GPX4', 'GADD45B'))
# DotPlot(pbmc.all.da, features=c('PDCLX', 'CD13', '' label = FALSE))
DotPlot(pbmc.all.da, features=c('SLC6A3', 'VAMP2', label = FALSE))
DotPlot(pbmc.all.da, features=c('SLC6A3', 'VAMP2', label = FALSE))

DotPlot(pbmc.all.da, features=c('TH', 'SLC6A3', label = FALSE))

table(pbmc.all.da$ASYN_range)



plot_aav_load(counts_t, "ASYN_BINS", "EGFR")

plot_aav_load(counts_t, "AAV-ASYN", "CHRNB3")

plot_aav_load(counts_t, "AAV-ASYN", "CHRNB3")
plot_aav_load(counts_t, "AAV-ASYN", "ATP2A3")
install.packages("pak") # if you have not installed "pak"
library("pak")
pak::pkg_install("pathfindR")
#############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")

devtools::install_github('cole-trapnell-lab/monocle3')

bfp.binned <- as.data.frame(pbmc.all.da$AAV_BFP)
head(bfp.binned)
# table(seurat.gfp.neg$AAV_BFP)
colnames(bfp.binned) <- c("BFP_count")
bfp.binned['BFP_range'] <- NA
bfp.binned = bfp.binned %>% mutate(BFP_range = case_when(
  BFP_count %in% 0 ~ "0",
  BFP_count %in% 1:10 ~ "1-3"))


bfp.binned = bfp.binned %>% mutate(BFP_range = case_when(
  BFP_count %in% 0 ~ "0",
  BFP_count %in% 1:93 ~ "1-93"))


table(bfp.binned$BFP_range)


# table(pbmc.all.da$N.ASYN.range)

# asyn.binned.df <- as.data.frame(asyn.binned)

pbmc.all.da <- AddMetaData(pbmc.all.da, bfp.binned["BFP_range"])
############################################################
Idents(object = pbmc.all.da) <- pbmc.all.da@meta.data$BFP_range
pbmc.all.da <- ScaleData(pbmc.all.da, scale = F)
DotPlot(pbmc.all.da, features=c('TH', 'SLC6A3', 'CALB1', 'CALB2'), scale = F)

table(Idents(pbmc.all.da))
#############################################################
#############################################################
pbmc.all.da.markers_0_1 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "1-3", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.da.markers_0_1, 40) -> top40
write.csv(pbmc.all.da.markers_0_1, paste0(dge_p, "markers_0_vs_1_bfp.csv"))
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40)))), scale = F) + coord_flip() + ggtitle("BFP 0 vs 1-93")
#############################################################
input_markers_bfp_bfp_0_vs_all <- pathfindr_df(pbmc.all.da.markers_0_1)
save_path = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/pathfindR_results"
output_markers_bfp_0_vs_1 <- run_pathfindR(input_markers_da_bfp_0_vs_all, gene_sets = "KEGG")

output_markers <- run_pathfindR(input_markers_da, gene_sets = "KEGG")


# output_markers_da_0_vs_93 <- run_pathfindR(input_markers_da_0_vs_93, output_dir = save_path, gene_sets = "KEGG")

term_gene_heatmap(output_markers_bfp_0_vs_1)
term_gene_graph(output_markers_bfp_0_vs_1)
table(Idents(pbmc.all.da))

pbmc.all.da.markers_2_93 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "2-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.da.markers_2_93, 40) -> top40
write.csv(pbmc.all.da.markers_2_93, paste0(dge_p, "markers_0_vs_2_93.csv"))
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()
#############################################################

pbmc.all.da.markers_1_2_93 <- FindMarkers(pbmc.all.da, ident.1 = "1", ident.2 = "2-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.da.markers_1_2_93, 40) -> top40
write.csv(pbmc.all.da.markers_1_2_93, paste0(dge_p, "markers_1_vs_2_93.csv"))
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()

#############################################################
#############################################################
# a-syn paracrine effect gfp- samples
pbmc.merged.gfp1 <- subset(x = pbmc.merged.final, subset = type == c("GFP-_asyn1"))
pbmc.merged.gfp2 <- subset(x = pbmc.merged.final, subset = type == c("GFP-_asyn2"))

seurat.gfp.neg <- merge(pbmc.merged.gfp1, pbmc.merged.gfp2)
table(seurat.gfp.neg$AAV_ASYN)

Idents(seurat.gfp.neg) <- seurat.gfp.neg$AAV_ASYN
asyn.binned <- as.data.frame(seurat.gfp.neg$AAV_ASYN)
head(asyn.binned)
colnames(asyn.binned) <- c("ASYN_count")
asyn.binned['ASYN_range'] <- NA
asyn.binned = asyn.binned %>% mutate(ASYN_range = case_when(
  ASYN_count %in% 0 ~ "0",
  ASYN_count %in% 1:20 ~ "1:20"))

head(asyn.binned)
table(asyn.binned$ASYN_range)
seurat.gfp.neg <- AddMetaData(seurat.gfp.neg, asyn.binned["ASYN_range"])

seurat.gfp.neg.markers <- FindAllMarkers(seurat.gfp.neg)

input_markers_gfp.neg <- pathfindr_df(seurat.gfp.neg.markers)

output_markers <- run_pathfindR(input_markers_gfp.neg, gene_sets = "KEGG")
