library(SingleCellExperiment)
library(Matrix)
library(scales)

setwd(".")

# How to read in 10X data for a single sample (output is a sparse matrix)

# Turn count matrix into a Seurat object (output is a Seurat object)
# The min.features argument specifies the minimum number of genes 
# that need to be detected per cell. This argument will filter out 
# poor quality cells that likely just have random barcodes encapsulated 
# without any cell present. Usually, cells with less than 100 genes detected 
# are not considered for analysis.
library(Seurat)

samples = c("Saline_Veh_A1","Saline_APAP_C1","Saline_3DDA_E1")
for (file in samples){
  seurat_data <- Read10X(data.dir = paste0(file, "/outs/filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, 
                                   min.features = 100, 
                                   project = file)
  assign(strsplit(file,"_")[[1]][3], seurat_obj)
}

rm(seurat_obj)
rm(seurat_data)
rm(file)

# Explore the metadata
# What do the columns of metadata mean?
#   
# orig.ident: this often contains the sample identity if known, but will default to "SeuratProject"
# nCount_RNA: number of UMIs per cell
# nFeature_RNA: number of genes detected per cell

# Reading in multiple samples
# Create each individual Seurat object for every sample

# Check the metadata in the new Seurat objects
# head(CFA@meta.data)
# head(Saline@meta.data)

# Create a merged Seurat object
merged_seurat <- merge(x = A1, 
                       y = list(C1,E1),
                       add.cell.id = c(strsplit(samples[1],"_")[[1]][2], 
                                       strsplit(samples[2],"_")[[1]][2], 
                                       strsplit(samples[3],"_")[[1]][2]))
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
# Explore merged metadata
View(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# Compute percent mito ratio
# merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
# merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
library(tidyverse)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$seq_folder, "^Saline_Veh_"))] <- "Veh"
metadata$sample[which(str_detect(metadata$seq_folder, "^Saline_APAP_"))] <- "APAP"
metadata$sample[which(str_detect(metadata$seq_folder, "^Saline_3DDA_"))] <- "3DDA"
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
result_dir = paste(strsplit(samples[1],"_")[[1]][3], "_", 
                   strsplit(samples[2],"_")[[1]][3], "_", 
                   strsplit(samples[3],"_")[[1]][3], sep="")
dir.create(result_dir, showWarnings = FALSE)
save(merged_seurat, file=paste(result_dir, "/merged_filtered_seurat.RData", sep=""))

# Visualize the number of cell counts per sample
# We see over 15,000 cells per sample, which is quite a bit more than the 12-13,000 expected. 
# It is clear that we likely have some junk 'cells' present.
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# Visualize the number UMIs/transcripts per cell
# The UMI counts per cell should generally be above 500, that is the low end of what we expect. 
# If UMI counts are between 500-1000 counts, 
# it is usable but the cells probably should have been sequenced more deeply.
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
# Visualize the distribution of genes detected per cell via histogram
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
# metadata %>%
#   ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
#   geom_point() +
#   scale_colour_gradient(low = "gray90", high = "black") +
#   stat_smooth(method=lm) +
#   scale_x_log10() +
#   scale_y_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 500) +
#   geom_hline(yintercept = 250) +
#   facet_wrap(~sample)
# Visualize the distribution of mitochondrial gene expression detected per cell
# metadata %>%
#   ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 0.2)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.85)


# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat,
                          subset= (nUMI >= 500) &
                            (nGene >= 250) &
                            (log10GenesPerUMI > 0.85))
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Create .RData object to load at any time
save(filtered_seurat, file=paste(result_dir, "/seurat_filtered.RData", sep=""))

# export_matrix <- as.matrix(filtered_seurat@assays$RNA@data)
# write.table(as.matrix(GetAssayData(object = filtered_seurat, slot = "counts")), 
#             'counts.txt', 
#             sep = '\t', row.names = T, col.names = T, quote = F)

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
# write.table(as.matrix(GetAssayData(object = seurat_phase, slot = "counts")), 
#             'counts_normal.csv', 
#             sep = ',', row.names = T, col.names = T, quote = F)
# Load cell cycle markers
# load("data/cycle.rda")

# Score cells for cell cycle
# seurat_phase <- CellCycleScoring(seurat_phase, 
                                 # g2m.features = g2m_genes, 
                                 # s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
# View(seurat_phase@meta.data)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
# DimPlot(seurat_phase,
#         reduction = "pca",
#         group.by= "Phase",
#         split.by = "Phase")

# adjust the limit for allowable object sizes within R (Default is 500 * 1024 ^ 2 = 500 Mb)
options(future.globals.maxSize = 4000 * 1024^2)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("Veh", "APAP", "3DDA")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  # split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=NULL, s.features=NULL)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = NULL)
}

# Check which assays are stored in objects
split_seurat$Veh@assays

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
library(RCurl)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
# Save integrated seurat object
saveRDS(seurat_integrated, paste(result_dir, "/integrated_seurat.rds", sep=""))
# seurat_integrated <- readRDS("A1_C1_E1/integrated_seurat.rds")
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")  
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)    
DimPlot(seurat_integrated,
        split.by = "sample")  

# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Explore resolutions
seurat_integrated@meta.data %>% 
  View()
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

## Calculation of UMAP
## DO NOT RUN (calculated in the last lesson)

# seurat_integrated <- RunUMAP(seurat_integrated, 
#                  reduction = "pca", 
#                  dims = 1:40)

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
# Explore whether clusters segregate by cell cycle phase
# DimPlot(seurat_integrated,
#         label = TRUE, 
#         split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
library(cowplot)
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("Nrxn3", "Lrrtm4"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("Kcnc2", "Lef1", "Ntng1","Synpo2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# FeaturePlot(seurat_integrated, 
#             reduction = "umap", 
#             features = c("MARCO", "ITGAM", "ADGRE1"), 
#             sort.cell = TRUE,
#             min.cutoff = 'q10', 
#             label = TRUE)
# FeaturePlot(seurat_integrated, 
#             reduction = "umap", 
#             features = c("FCER1A", "CST3"), 
#             sort.cell = TRUE,
#             min.cutoff = 'q10', 
#             label = TRUE)
# FeaturePlot(seurat_integrated, 
#             reduction = "umap", 
#             features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
#             sort.cell = TRUE,
#             min.cutoff = 'q10', 
#             label = TRUE)

## DO NOT RUN THIS CODE ##

# Find markers for every cluster compared to all remaining cells, report only the positive ones
# markers <- FindAllMarkers(object = seurat_integrated, 
#                           only.pos = TRUE,
#                           logfc.threshold = 0.25)

DefaultAssay(seurat_integrated) <- "RNA"
# install.packages('BiocManager')
# BiocManager::install('multtest')
# install.packages('metap')
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    # left_join(y = unique(annotations[, c("gene_name", "description")]),
    #           by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(1,2), get_conserved)

# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (Veh_avg_log2FC + APAP_avg_log2FC + `3DDA_avg_log2FC`) /3) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
View(top10)

# Plot interesting marker gene expression for cluster 1
FeaturePlot(object = seurat_integrated, 
            features = c("Klhl1", "Ebf1", "Cacna2d1", "Cntnap5b", "Cntnap5c"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

# Vln plot - cluster 1
VlnPlot(object = seurat_integrated, 
        features = c("Klhl1", "Ebf1", "Cacna2d1", "Cntnap5b", "Cntnap5c"))

# Determine differentiating markers for CD4+ T cell
cd4_tcells <- FindMarkers(seurat_integrated,
                          ident.1 = 2,
                          ident.2 = c(0,4,10,18))                  

# Add gene symbols to the DE table
# cd4_tcells <- cd4_tcells %>%
#   rownames_to_column(var = "gene") %>%
#   left_join(y = unique(annotations[, c("gene_name", "description")]),
#             by = c("gene" = "gene_name"))

# Reorder columns and sort by padj      
cd4_tcells <- cd4_tcells[, c(1, 3:5,2)]

cd4_tcells <- cd4_tcells %>%
  dplyr::arrange(p_val_adj) 

# View data
View(cd4_tcells)

# Rename all identities
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "Naive or memory CD4+ T cells",
                                  "1" = "CD14+ monocytes",
                                  "2" = "Naive or memory CD4+ T cells",
                                  "3" = "CD14+ monocytes",
                                  "4" = "CD4+ T cells",
                                  "5" = "CD8+ T cells",
                                  "6" = "B cells",
                                  "7" = "Stressed cells / Activated T cells",
                                  "8" = "NK cells",
                                  "9" = "FCGR3A+ monocytes",
                                  "10" = "CD4+ T cells",
                                  "11" = "B cells",
                                  "12" = "NK cells",
                                  "13" = "CD8+ T cells",
                                  "14" = "CD14+ monocytes",
                                  "15" = "Conventional dendritic cells",
                                  "16" = "Megakaryocytes",
                                  "17" = "B cells", 
                                  "18" = "CD4+ T cells", 
                                  "19" = "Plasmacytoid dendritic cells", 
                                  "20" = "Mast cells")

# Plot the UMAP
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
