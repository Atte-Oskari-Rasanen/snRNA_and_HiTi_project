library(Seurat)
library(Matrix)

##############################
#Ref own
p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/anndata_from_loom"
# p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref/out"

mat <- read.csv(paste(p,"COUNT_MATRIX_ADATA_1000.csv", sep="/"), row.names=NULL)
# mat <- read.csv(paste(p,"ad_unfilt_matrix.csv", sep="/"), row.names=NULL)

cells <- read.csv(paste(p,"CELLS_ADATA_1000.csv", sep="/"), row.names=NULL)
# cells <- read.csv(paste(p,"/ad_unfilt_cells.csv", sep="/"), row.names=NULL)

genes <- read.csv(paste(p,"GENES_ADATA_1000.csv", sep="/"), row.names=NULL)
# genes <- read.csv(paste(p,"ad_unfilt_genes.csv", sep="/"), row.names=NULL)

cellids <- read.csv(paste(p, "CELL_IDs_ADATA_1000.csv", sep="/"), row.names=NULL)
cells2 <- unlist(cells[2])
cellids2 <- unlist(cellids[2])
genes2 <- unlist(genes[2])
cells2 <- head(cells2, -6)
cellids2 <- head(cellids2, -6)

names(cells2) <- gsub("[[:digit:]]", "", names(cells2))
names(cells2) <- gsub("[[:digit:]]", "", names(cellids2))
names(genes2) <- gsub("[[:digit:]]", "", names(genes2))
dim(mat)
mat <-mat[-27999,-11011]
mat <- t(mat)
########################################
#rownames(mat) <- genes2
# colnames(mat) <- cellids2
mat <- as(as.matrix(mat), "sparseMatrix")
libsizes <- colSums(mat)
size.factors <- libsizes/mean(libsizes)
logcounts <- log2(t(t(mat)/size.factors) + 1)
logcounts[1:3,1:3]
se <- SummarizedExperiment(list(logcounts=logcounts))
rownames(se) <- genes2
colnames(se) <- cellids2
se$label.main <- cells2
library("SingleCellExperiment")

#################################################
ann_p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_df.csv"
scifi5.df <- read.csv(ann_p)

genes <- scifi5.df[,1]
genes
scifi5.df_norows <- scifi5.df[,-1]
scifi5.df[1:3,1:3]
scifi5.matrix <- as.matrix(scifi5.df_norows)
sparse.scifi5.matrix <- as(scifi5.matrix, "sparseMatrix")
rownames(sparse.scifi5.matrix) <- genes

scifi5.matrix[1:3,1:3]

seurat.scifi <- CreateSeuratObject(counts = sparse.scifi5.matrix, min.cells = 3, min.features  = 100, project = "Scifi5", assay = "RNA")

df <- t(scifi5.df)
colnames(df) <- genes
df[1:3,1:3]
tail(colnames(df))
df[c("AAV_ASYN")]
scifi5.df[["AAV_ASYN"]]

colnames(df)

aav_meta <- df[,c("AAV_ASYN", "AAV_BFP")]
aav_meta1 <- df[,c("AAV_ASYN")]
aav_meta2 <- df[,c("AAV_BFP")]

aav_meta[is.na(aav_meta)] = 0
head(aav_meta)
head(scifi5.df)

scifi5.df[1:3,1:3]



seurat.scifi@meta.data$AAV <- aav_meta
seurat.scifi <- AddMetaData(seurat.scifi, aav_meta1, col.name = "AAV_ASYN")
seurat.scifi <- AddMetaData(seurat.scifi, aav_meta2, col.name = "AAV_BFP")

# scifi.data <- Read10X(data.dir = "/media/data/AtteR/projects/scifi-analysis/outputs_starsolo/Scifi5/Scifi_library_2_S2_SciFi5_oDT/GeneFull/filtered")

head(scifi.data)
sc

# seurat.scifi2 <- CreateSeuratObject(counts = scifi.data, min.cells = 3, min.features  = 200, project = "Scifi5", assay = "RNA")

VlnPlot(seurat.scifi, features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI"), ncol = 2)

FeatureScatter(seurat.scifi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seurat.scifi <- subset(seurat.scifi, subset = nFeature_RNA > 200 & nFeature_RNA < 600 & nCount_RNA < 26000)

seurat.scifi <- NormalizeData(seurat.scifi, normalization.method = "LogNormalize", scale.factor = 1000)

seurat.scifi <- FindVariableFeatures(seurat.scifi, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.scifi), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.scifi)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2

###########################################################
###########################################################
seurat.scifi@assays$RNA

counts_all <- as.data.frame(seurat.scifi$nCount_RNA)
library (plyr)
df_c <- ldply (counts_all, data.frame)
df1 <- df_c[1:length(rownames(df_c)),2]
#df1
p1 <- hist(df1, main="", 
           , xlab = "Number of UMIs per cell Scifi5",col = "gray",border = "black", breaks=100)
counts_all <- as.data.frame(seurat.scifi$nFeature_RNA)

df_c <- ldply (counts_all, data.frame)
df1 <- df_c[1:length(rownames(df_c)),2]
#df1
hist(df1, main="", 
     , xlab = "Number of genes per cell Scifi5",col = "gray",border = "black", breaks=100)
###########################################################
###########################################################

all.genes <- rownames(seurat.scifi)
seurat.scifi <- ScaleData(seurat.scifi, features = all.genes)

seurat.scifi <- RunPCA(seurat.scifi, features = VariableFeatures(object = seurat.scifi))
ElbowPlot(seurat.scifi)
seurat.scifi <- FindNeighbors(seurat.scifi, dims = 1:30)
seurat.scifi <- FindClusters(seurat.scifi, resolution = 1.2)

seurat.scifi <- RunUMAP(seurat.scifi, dims = 1:30)

DimPlot(seurat.scifi, reduction = "umap")
head(seurat.scifi[[]])

counts_s4 <- seurat.scifi$RNA@counts
counts <- as.data.frame(counts_s4)
counts[1:5,1:5]
write.csv(counts, "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/da_subtype_counts_scifi.csv")

final_ann <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/DA_subtypes_scifi_meta.csv"

anns_with_da_subtypes_final <- read.csv(final_ann)
head(anns_with_da_subtypes_final)
colnames(anns_with_da_subtypes_final) <- "cell_types"
table(anns_with_da_subtypes_final[,2])

anns_with_da_subtypes_final[1:3,2]

anns_with_da_subtypes_final[1:3,1:2]
anns_with_da_subtypes_final2 <- as.data.frame(anns_with_da_subtypes_final[,2])
colnames(anns_with_da_subtypes_final2) <- "cell_types"
rownames(anns_with_da_subtypes_final2) <- anns_with_da_subtypes_final[,1]


seurat.scifi <- AddMetaData(seurat.scifi, anns_with_da_subtypes_final2, col.name = "Annotations_specific")



Idents(seurat.scifi) <- seurat.scifi@meta.data$Annotations_specific
seurat.scifi$AAV_ASYN
library("ggplot2")
FeaturePlot(seurat.scifi, reduction = "umap", "AAV_ASYN", label = FALSE) + ggtitle("AAV_ASYN copy number")

p3 <- FeaturePlot(seurat.scifi, reduction = "umap", "Th", label = FALSE)
p3

asyn.binned <- as.data.frame(seurat.scifi$AAV_ASYN)
colnames(asyn.binned) <- "ASYN_count"
head(asyn.binned)

asyn.binned["ASYN_count"] <- lapply(asyn.binned["ASYN_count"], strtoi)
asyn.binned['ASYN_range'] <- NA
asyn.binned = asyn.binned %>% mutate(ASYN_range = case_when(
  ASYN_count %in% 0 ~ "0",
  ASYN_count %in% 1:5 ~ "1:5",
  ASYN_count %in% 6:14 ~ "6:14"))


table(asyn.binned$ASYN_range)
max(asyn.binned$ASYN)

seurat.scifi.markers_0_1_14 <- FindMarkers(seurat.scifi, ident.1 = "0", ident.2 = "1:14", logfc.threshold = 0.15, test.use = "wilcox", min.pct = 0.1)

#############################################################
seurat.scifi <- AddMetaData(seurat.scifi, asyn.binned)
table(seurat.scifi$ASYN_range)
seurat.scifi@active.ident <- factor(seurat.scifi@active.ident, 
                            levels=c("0", "1", "2:14"))
############################################################
Idents(seurat.scifi) <- seurat.scifi$ASYN_range
table(Idents(seurat.scifi))

pbmc.all.da.scifi <- FindAllMarkers(seurat.scifi)
pbmc.all.da.scifi

bpf.binned <- seurat.scifi$AAV_BFP
bpf.binned[bpf.binned >= 1 & bpf.binned <= 5] <- "1-5"
bpf.binned[bpf.binned > 5 & bpf.binned < 11] <- "6-10"
bpf.binned[bpf.binned > 10 & bpf.binned <= 20] <- "11-20"
bpf.binned[bpf.binned > 20 & bpf.binned <= 30] <- "21-30"
bpf.binned[bpf.binned > 30 & bpf.binned <= 40] <- "31-40"
bpf.binned[bpf.binned > 40 & bpf.binned <= 93] <- "41-93"
############################################################
bpf.binned.df <- as.data.frame(bpf.binned)
colnames(bpf.binned.df) <- c("N.BPF.range")
rownames(bpf.binned.df) <- rownames(seurat.scifi[[]])
head(bpf.binned.df)
unique(bpf.binned.df)
seurat.scifi <- AddMetaData(seurat.scifi, bpf.binned.df)

table(seurat.scifi$N.BPF.range)

table(Idents(pbmc.all.da))

seurat.scifi[["old.ident"]] <- Idents(object = seurat.scifi)


rds_out <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_proc.rds"

saveRDS(seurat.scifi, rds_out)
seurat.scifi <- readRDS(rds_out)

DimPlot(seurat.scifi)
Idents(seurat.scifi) <- seurat.scifi$ASYN_range
DotPlot(seurat.scifi, features=c('Th', 'Slc6a3', 'Calb1', 'Calb2'), scale = F)


seurat.scifi <- readRDS(rds_out)
hist(pbmc.all.da$AAV_ASYN, main="Injected Asyn copy number in cells",
     xlab="Copy number of Asyn",
     xlim=c(0,50), breaks=50)

Idents(object = seurat.scifi) <- seurat.scifi@meta.data$N.ASYN.range
table(Idents(seurat.scifi))
head(pbmc.all.da.markers_0_1.5, 10)
seurat.scifi.markers <- FindAllMarkers(seurat.scifi)

# could not find features past logfc threshold

seurat.scifi.markers_0_1 <- FindMarkers(seurat.scifi, ident.1 = "0", ident.2 = "1", logfc.threshold = 0.15, test.use = "wilcox", min.pct = 0.1)

head(seurat.scifi.markers_0_1, 40) -> top40
# seurat.scifi.markers_0_1_5 %>% top_n(-40, p_val_adj) -> top40
pbmc.scaled <- ScaleData(seurat.scifi, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()
#############################################
head(seurat.scifi.markers_0_1, 10)
seurat.scifi.markers_0_2_14 <- FindMarkers(seurat.scifi, ident.1 = "0", ident.2 = "2:14")
head(pbmc.all.da.markers_1.5_11.20, 10)

head(pbmc.all.da.markers_1.5_11.20, 40) -> top40
# seurat.scifi.markers_0_1_5 %>% top_n(-40, p_val_adj) -> top40
pbmc.scaled <- ScaleData(seurat.scifi, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()
#############################################
#############################################
head(pbmc.all.da.markers_1.5_10, 10)
pbmc.all.da.markers_1.5_11.20 <- FindMarkers(seurat.scifi, ident.1 = "1-5", ident.2 = "11-20")
head(pbmc.all.da.markers_1.5_11.20, 10)

head(pbmc.all.da.markers_1.5_11.20, 40) -> top40
# seurat.scifi.markers_0_1_5 %>% top_n(-40, p_val_adj) -> top40
pbmc.scaled <- ScaleData(seurat.scifi, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()

###################################################
Idents(object = seurat.scifi) <- seurat.scifi@meta.data$N.BPF.range

# No significant difference found between groups in terms of gene expression

pbmc.all.da.markers_0_1.5_bfp <- FindMarkers(seurat.scifi, ident.1 = "0", ident.2 = "1-5", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.2)
head(pbmc.all.da.markers_0_1.5_bfp, 40) -> top40
write.csv(pbmc.all.da.markers_0_1.5_bfp, paste0(dge_p, "markers_0_vs_15_bfp.csv"))
pbmc.scaled <- ScaleData(seurat.scifi, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40))))
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip()


pbmc.all.da.copy@meta.data$N.BPF.range
#############################################################
# subset into own data set based on asyn range. then find highly variable features
seurat.scifi_0asyn <- subset(x = seurat.scifi, subset = AAV_ASYN == c("0"))
seurat.scifi_0asyn <- FindVariableFeatures(seurat.scifi_0asyn, selection.method = "vst", nfeatures = 2000)


p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_var_genes/"
# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(seurat.scifi_0asyn), 50)
top50
write.csv(top50, paste0(p,"scifi5_asyn0.csv"))
top10
#############################################################
seurat.scifi_1_5asyn <- subset(x = seurat.scifi, subset = AAV_ASYN == c("1-5"))
seurat.scifi_1_5asyn <- FindVariableFeatures(seurat.scifi_1_5asyn, selection.method = "vst", nfeatures = 2000)
p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_var_genes/"
# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(seurat.scifi_1_5asyn), 50)
top50
write.csv(top50, paste0(p,"scifi5_asyn1_5.csv"))
top10

seurat.scifi_1_5asyn <- subset(x = seurat.scifi, subset = AAV_ASYN == c("0"))
seurat.scifi_0asyn <- FindVariableFeatures(seurat.scifi_0asyn, selection.method = "vst", nfeatures = 2000)
#############################################################
#############################################################
seurat.scifi_0_1_5bfp <- subset(x = seurat.scifi, subset = AAV_ASYN == c("1-5"))
seurat.scifi_1_5asyn <- FindVariableFeatures(seurat.scifi_1_5asyn, selection.method = "vst", nfeatures = 2000)
p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_var_genes/"
# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(seurat.scifi_1_5asyn), 50)
top50
write.csv(top50, paste0(p,"scifi5_asyn1_5.csv"))
top10

seurat.scifi_1_5asyn <- subset(x = seurat.scifi, subset = AAV_ASYN == c("0"))
seurat.scifi_0asyn <- FindVariableFeatures(seurat.scifi_0asyn, selection.method = "vst", nfeatures = 2000)
#############################################################

p <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_var_genes/"
# Identify the 10 most highly variable genes
top50 <- head(VariableFeatures(seurat.scifi_0asyn), 50)
top50
write.csv(top50, paste0(p,"scifi5_asyn0.csv"))
top10



sce.scifi <- as.SingleCellExperiment(seurat.scifi)
library("SingleR")

# would need to subset based on DA markers but we did not do this 
# with the parse stuff so what justifies us to do this now?

pred <- SingleR(test=sce.scifi, ref=se, labels=se$label.main)

table(pred$labels)
seurat.scifi$annotations <- pred$labels
p = "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/annotations_final_nonspecific_da_scifi.csv"

write.csv(seurat.scifi$annotations, p)

          
Idents(object = seurat.scifi) <- seurat.scifi$annotations

seurat.scifi.da <- subset(x = seurat.scifi, subset = annotations == c("Dopaminergic neurons, ventral midbrain (SNc, VTA)"))
counts_s4 <- seurat.scifi.da$RNA@counts
counts <- as.data.frame(counts_s4)
counts[1:5,1:5]
write.csv(counts, "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/da_subtype_counts_scifi.csv")

rds_out <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_ann_preproc.rds"
rds_out <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_da_sub_ann_preproc.rds"
rds_out <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/scifi5_proc.rds"

saveRDS(seurat.scifi, rds_out)
seurat.scifi <- readRDS(rds_out)
###########################################
table(seurat.scifi$Annotations_specific)

final_ann <- "/media/data/AtteR/projects/parsebio/ParsePipeline/docker_singlecell/analysis_stuff/ref_files_final/final_annotations_scifi_with_dasubtypes.csv"
anns_with_da_subtypes_final <- read.csv(final_ann, row.names="genes")
head(anns_with_da_subtypes_final)


seurat.scifi$Annotations_specific <- anns_with_da_subtypes_final
seurat.scifi <- AddMetaData(seurat.scifi, anns_with_da_subtypes_final, col.name = "Annotations_specific")

table(seurat.scifi$Annotations_specific)
Idents(seurat.scifi) <- seurat.scifi@meta.data$Annotations_specific
table(Idents(seurat.scifi))
DimPlot(seurat.scifi)


library(ggplot2)
table(seurat.scifi$ASYN_range)

seurat.scifi

Idents(seurat.scifi) <- seurat.scifi$ASYN_range
Idents(seurat.scifi) <- seurat.scifi@meta.data$seurat_clusters

Idents(seurat.scifi) <- seurat.scifi$ASYN_range
table(Idents(seurat.scifi))

pbmc.all.da.scifi <- FindAllMarkers(seurat.scifi)
pbmc.all.da.scifi

pbmc.all.da.markers_0_1 <- FindMarkers(seurat.scifi, ident.1 = "0", ident.2 = "2", logfc.threshold = 0.1, test.use = "wilcox", min.pct = 0.1)

head(pbmc.all.da.markers_0_1, 10)
write.csv(pbmc.all.da.markers_0_1, paste0(dge_p, "markers_0_vs_1.csv", ))

head(pbmc.all.da.markers_0_1, 40) -> top40
pbmc.scaled <- ScaleData(pbmc.all.da, features = as.character(unique(rownames(top40))))
DoHeatmap(pbmc.scaled, features = as.character(unique(rownames(top40)))) + ggtitle("ASyn 0 vs 1")
DotPlot(pbmc.scaled, features = rev(as.character(unique(rownames(top40))))) + coord_flip() + ggtitle("ASyn 0 vs 1")
