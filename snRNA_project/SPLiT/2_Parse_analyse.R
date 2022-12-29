library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
# install.packages(Seurat)
# BiocManager::install("slingshot")
# library(slingshot)

out_seurat <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/seurat_objects/"
# saveRDS(pbmc.merged.final,paste0(out_seurat, "Seurat_Parse_Full_Ann.rds"))
pbmc.merged.final <- readRDS(paste0(out_seurat, "Seurat_Parse_Full_Ann.rds"))
table(pbmc.merged.final$Annotations_specific)
p1 <- DimPlot(pbmc.merged.final)
unique(pbmc.merged.final$AAV_ASYN)
p2 <- FeaturePlot(pbmc.merged.final, pbmc.merged.final$AAV_ASYN, label = FALSE) + ggtitle("AAV_ASYN copy number")
p1 | p2

#############################################
# Histograms showing the number of different cell types

pbmc.merged.gfp1 <- subset(x = pbmc.merged.final, subset = type == c("GFP-_asyn1"))
pbmc.merged.gfp2 <- subset(x = pbmc.merged.final, subset = type == c("GFP-_asyn2"))
seurat.gfp.neg <- merge(pbmc.merged.gfp1, y=pbmc.merged.gfp2)

seurat.gfp.pos.asyn <- subset(x = pbmc.merged.final, subset = type == c("GFP+_asyn"))

seurat.gfp.pos.bfp <- subset(x = pbmc.merged.final, subset = type == c("GFP+_bfp"))
anns_gfp_neg <- as.data.frame(table(seurat.gfp.neg$Annotations_specific))
head(anns)
colnames(anns) <- c('Cell_type','Count')


initial_anns <- as.data.frame(table(pbmc.merged.final$annotations_1000))
colnames(initial_anns) <- c('Cell_type','Count')

final_anns <- as.data.frame(table(pbmc.merged.final$Annotations_specific))
colnames(final_anns) <- c('Cell_type','Count')


anns_neg <- as.data.frame(table(seurat.gfp.neg$Annotations_specific))
head(anns_neg)
colnames(anns_neg) <- c('Cell_type','Count')

anns_asyn_pos <- as.data.frame(table(seurat.gfp.pos.asyn$Annotations_specific))
head(anns_asyn_pos)
colnames(anns_asyn_pos) <- c('Cell_type','Count')

anns_bfp_pos <- as.data.frame(table(seurat.gfp.pos.bfp$Annotations_specific))
head(anns)
colnames(anns_bfp_pos) <- c('Cell_type','Count')

ggplot(anns_neg, aes(Cell_type, Count)) +
  geom_col() + 
  coord_flip() + ylab("Number of cells") + xlab("Cell types GFP- a-syn") + theme(text = element_text(size = 20)) 



# + geom_point
  # scale_x_continuous(br() +eaks = scales::pretty_breaks(n = 10))
  # scale_x_continuous(limits = c(0, 1000))
ggplot(anns_asyn_pos, aes(Cell_type, Count)) +
  geom_col() + 
  coord_flip() + ylab("Number of cells") + xlab("Cell types GFP+ a-syn") + theme(text = element_text(size = 20)) 

ggplot(anns_bfp_pos, aes(Cell_type, Count)) +
  geom_col() + 
  coord_flip() + ylab("Number of cells") + xlab("Cell types GFP+ tagBFP") + theme(text = element_text(size = 20)) 


ggplot(initial_anns, aes(Cell_type, Count)) +
  geom_col() + 
  coord_flip() + ylab("Number of cells") + xlab("Cell types all") + theme(text = element_text(size = 20)) 

ggplot(initial_anns, aes(Cell_type, Count)) +
  geom_col() + 
  coord_flip() + ylab("Number of cells") + xlab("Cell types all SingleR") + theme(text = element_text(size = 20)) 


ggplot(final_anns, aes(Cell_type, Count)) +
  geom_col() + 
  coord_flip() + ylab("Number of cells") + xlab("Cell types all final annotations") + theme(text = element_text(size = 20)) 

#############################################
# Pie chart
slices <- c(10, 12, 4, 16, 8)
lbls <- as.vector(seurat.gfp.neg$Annotations_specific)
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries") 

table(pbmc.merged.final@meta.data$type)

anns <- as.data.frame(table(seurat.gfp.neg$Annotations_specific))
head(anns)
colnames(anns) <- c('Cell_type','Percentage')


library(ggplot2)
ggplot(DF, aes(x = Condition, y = percent, fill = Cell_Type))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste(percent,"%")), position = position_stack(vjust = 0.5))

library(ggrepel)
#######################################
# Dotplots for different markers across cell types
Idents(object = pbmc.merged.final) <- pbmc.merged.final@meta.data$Annotations_specific
DotPlot(pbmc.merged.final, features=c('GAD1', 'GAD2', 'PAX2', 'SLC32A1', label = FALSE))
DotPlot(pbmc.merged.final, features=c('SLC17A7', 'SLC17A6', 'MEIS2','SLC1A1', label = FALSE))
DotPlot(pbmc.merged.final, group.by = "final_annotations", features=c('COL6A2', 'VIM', 'PDGFRB', label = FALSE))

######################################################
pbmc.merged.final
cells.use <- colnames(pbmc.merged.final)[which(pbmc.merged.final[[]]['sample'] == "619_GFP")]
pbmc.merged_bfp <- subset(pbmc.merged.final, cells = cells.use)
table(pbmc.merged_bfp$final_annotations)
######################################################
cells.use <- colnames(pbmc.merged.final)[which(pbmc.merged.final[[]]['sample'] == "650_GFP")]
pbmc.merged_asyn <- subset(pbmc.merged.final, cells = cells.use)
table(pbmc.merged_asyn$final_annotations)
######################################################
cells.use <- colnames(pbmc.merged.final)[which(pbmc.merged.final[[]]['sample'] == "637_noGFP")]
pbmc.merged_nogfp1 <- subset(pbmc.merged.final, cells = cells.use)
table(pbmc.merged_nogfp1$final_annotations)
######################################################
cells.use <- colnames(pbmc.merged.final)[which(pbmc.merged.final[[]]['sample'] == "648_noGFP")]
pbmc.merged_nogfp2 <- subset(pbmc.merged.final, cells = cells.use)
table(pbmc.merged_nogfp2$final_annotations)

pbmc.all.da <- merge(pbmc.merged_bfp, pbmc.merged_asyn)
table(pbmc.all.da$Annotations_specific)
pbmc.all.nonda <- merge(pbmc.merged_nogfp1, pbmc.merged_nogfp2)
table(pbmc.all.nonda$Annotations_specific)

# take the non da cells from the whole seurat object and continue with their analysis
pbmc.merged.nonda <- pbmc.merged.final[, !(pbmc.merged.final$Annotations_specific %in% c("DA_SNc", "DA_VTA1", "DA_VTA3", "DA_VTA4", "Dopaminergic neurons, ventral midbrain (SNc, VTA)"))]
table(pbmc.merged.nonda$Annotations_specific)

##########################################################
# Bin into a-syn groups and set as identities for the cells, then perform DGE analysis, visualise with dotplots and heatmaps
asyn.binned <- as.data.frame(pbmc.merged.nonda$AAV_ASYN)
head(asyn.binned)
colnames(asyn.binned) <- c("ASYN_count")
asyn.binned['ASYN_range'] <- NA
asyn.binned = asyn.binned %>% mutate(ASYN_range = case_when(
  ASYN_count %in% 0 ~ "0",
  ASYN_count %in% 1 ~ "1",
  ASYN_count %in% 2:3 ~ "2-3",
  ASYN_count %in% 4:93 ~ "4-93"))

head(asyn.binned)
table(asyn.binned$ASYN_range)

# table(pbmc.all.da$N.ASYN.range)

# asyn.binned.df <- as.data.frame(asyn.binned)
rownames(asyn.binned) <- rownames(pbmc.merged.nonda)
head(asyn.binned)
unique(asyn.binned)

pbmc.merged.nonda <- AddMetaData(pbmc.merged.nonda, asyn.binned["ASYN_range"])
pbmc.merged.nonda$ASYN_range
table(pbmc.merged.nonda$ASYN_range)

Idents(pbmc.merged.nonda) <- pbmc.merged.nonda$ASYN_range
DotPlot(pbmc.merged.nonda, features=c('TH', 'SLC6A3', 'CALB1', 'CALB2'), scale = F)
pbmc.all.da.markers_0_1 <- FindMarkers(pbmc.merged.nonda, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.nonda.markers_0_1, 10)
write.csv(pbmc.all.nonda.markers_0_1, paste0(dge_p, "markers_0_vs_1.csv", ))

head(pbmc.all.nonda.markers_0_1, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.merged.nonda, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 0 vs 1")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 0 vs 1")


pbmc.all.nonda.markers_0_93 <- FindMarkers(pbmc.merged.nonda, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.nonda.markers_0_93, 10)
write.csv(pbmc.all.nonda.markers_0_93, paste0(dge_p, "markers_0_vs_1.csv", ))

head(pbmc.all.nonda.markers_0_93, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.merged.nonda, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 0 vs 93")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 0 vs 93")

#################################################################
pbmc.all.nonda.markers_1_93 <- FindMarkers(pbmc.merged.nonda, ident.1 = "1", ident.2 = "4-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.nonda.markers_1_93, 10)
write.csv(pbmc.all.nonda.markers_1_93, paste0(dge_p, "markers_0_vs_1.csv", ))

head(pbmc.all.nonda.markers_1_93, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.merged.nonda, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 1 vs 4-93")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 1 vs 4-93")

# Further subset into DA types
pbmc.all.da.snc <- subset(x = pbmc.all.da, subset = Annotations_specific2 == c("DA_SNc"))
pbmc.all.da.vta1 <- subset(x = pbmc.all.da, subset = Annotations_specific2 == c("DA_VTA1"))
pbmc.all.da.vta3 <- subset(x = pbmc.all.da, subset = Annotations_specific2 == c("DA_VTA3"))
pbmc.all.da.vta4 <- subset(x = pbmc.all.da, subset = Annotations_specific2 == c("DA_VTA4"))
pbmc.all.da.all <- subset(x = pbmc.all.da, subset = Annotations_specific2 == c("Dopaminergic neurons, ventral midbrain (SNc, VTA)"))
pbmc.all.da <- merge(x = pbmc.all.da.snc, y = c(pbmc.all.da.vta1, pbmc.all.da.vta3, pbmc.all.da.all))

table(pbmc.all.da$Annotations_specific)

##########################################################
asyn.binned <- as.data.frame(pbmc.all.da$AAV_ASYN)
head(asyn.binned)
colnames(asyn.binned) <- c("ASYN_count")
asyn.binned['ASYN_range'] <- NA
asyn.binned = asyn.binned %>% mutate(ASYN_range = case_when(
  ASYN_count %in% 0 ~ "0",
  ASYN_count %in% 1 ~ "1",
  ASYN_count %in% 2:3 ~ "2-3",
  ASYN_count %in% 4:8 ~ "4-8",
  ASYN_count %in% 9:15 ~ "9-15", ASYN_count %in% 16:93 ~ "16-93"))

head(asyn.binned)
table(asyn.binned$ASYN_range)

# table(pbmc.all.da$N.ASYN.range)

# asyn.binned.df <- as.data.frame(asyn.binned)
rownames(asyn.binned) <- rownames(pbmc.all.da)
head(asyn.binned)
unique(asyn.binned)

pbmc.all.da <- AddMetaData(pbmc.all.da, asyn.binned["ASYN_range"])
pbmc.all.da$ASYN_range
table(pbmc.all.da$ASYN_range)

pbmc.all.da$meta.data
pbmc.all.da <- FindVariableFeatures(pbmc.all.da)
pbmc.all.da <- ScaleData(pbmc.all.da)
pbmc.all.da <- RunPCA(pbmc.all.da)
pbmc.all.da <- RunUMAP(pbmc.all.da, dims=1:30, reduction="pca")
pbmc.all.da <- FindNeighbors(pbmc.all.da, reduction = "pca", dims = 1:20) #10
pbmc.all.da <- FindClusters(pbmc.all.da, resolution = 0.5)
DimPlot(pbmc.all.da, group.by = "ASYN_range")


cluster_data = pbmc.all.da@meta.data$seurat_clusters
da_count_mtx <- as.data.frame(pbmc.all.da@assays$originalexp@counts)
da_count_mtx[1:3,1:3]

save_p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/aav_groupings_data_viaenv/"
write.csv(da_count_mtx, paste0(save_p, "da_count_mtx.csv"))
write.csv(asyn.binned["ASYN_range"], paste0(save_p, "aav_asyn_da.csv"))
write.csv(cluster_data, paste0(save_p, "seurat_clusters.csv"))

Idents(pbmc.all.da) <- pbmc.all.da$ASYN_range
DotPlot(pbmc.all.da, features=c('TH', 'SLC6A3', 'CALB1', 'CALB2'), scale = F)
pbmc.all.da <- RunPCA(pbmc.all.da)
pbmc.all.da <- RunUMAP(pbmc.all.da, dims = 1:30)
##########################################################

pbmc.all.da[["old.ident"]] <- Idents(object = pbmc.all.da)

hist(pbmc.all.da$AAV_ASYN, main="Injected Asyn copy number in cells",
     xlab="Copy number of Asyn",
     xlim=c(0,50), breaks=50)
outp
saveRDS(pbmc.all.da, "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/pbmc_parse_da_subset.rds")



##################################################

pbmc.all.da <- readRDS("/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/pbmc_parse_da_subset.rds")

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
head(pbmc.all.da.markers_0__9_93)
log_change <- as.vector(pbmc.all.da.markers_0__9_93$avg_log2FC)
p_val <- as.vector(pbmc.all.da.markers_0__9_93$p_val)
gene_names <- as.vector(pbmc.all.da.markers_0__9_93)

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
# Plot the influence AAV load on a selected gene
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


plot_aav_load(counts_t, "AAV-ASYN", "CHRNB3")
plot_aav_load(counts_t, "AAV-ASYN", "CHRNB3")
plot_aav_load(counts_t, "AAV-ASYN", "ATP2A3")
install.packages("pak") # if you have not installed "pak"
library("pak")
pak::pkg_install("pathfindR")
#############################################################
# Repeat the procedure of binning and DGE with tagBFP control
bfp.binned <- as.data.frame(pbmc.all.da$AAV_BFP)
head(bfp.binned)
colnames(bfp.binned) <- c("BFP_count")
bfp.binned['BFP_range'] <- NA

bfp.binned = bfp.binned %>% mutate(BFP_range = case_when(
  BFP_count %in% 0 ~ "0",
  BFP_count %in% 1 ~ "1",
  BFP_count %in% 2:93 ~ "2-93"))

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
pbmc.all.da.markers_0_1 <- FindMarkers(pbmc.all.da, ident.1 = "0", ident.2 = "1-93", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.3)
head(pbmc.all.da.markers_0_1, 40) -> top40
write.csv(pbmc.all.da.markers_0_1, paste0(dge_p, "markers_0_vs_1_bfp.csv"))
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40)))), scale = F) + coord_flip() + ggtitle("BFP 0 vs 1-93")
#############################################################
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
