# Title: 05 - Subsetting Epithelial Cells
# Date: 2021Sept28
# Author: Jayne Wiarda

library(Seurat)
library(SeuratDisk)
library(ggplot2)

gut <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpAll.h5Seurat')
DefaultAssay(gut) <- 'integrated'

# Re-identify PCdims:
pct <- gut[["pca"]]@stdev / sum(gut[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
#[1] 15
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps

# Perform cell clustering
gut <- FindClusters(gut, dims = PCdims, reduction = 'pca') # We already calculated FindNeighbors previously, so just run FindClusters & with default resolution
DimPlot(gut, group.by = 'seurat_clusters', label = TRUE) # look at clusters on UMAP
DimPlot(gut, group.by = 'seurat_clusters', label = FALSE, split.by = 'orig.ident', ncol = 2) & NoAxes() & NoLegend() # look at clusters on UMAP split by sample IDs

# Perform hierarchical clustering of cell clusters
DefaultAssay(gut) <- 'RNA'
gut <- BuildClusterTree(gut, 
                       dims = PCdims, 
                       assay = "PCA")
data.tree <- Tool(object =gut, 
                  slot = "BuildClusterTree") 
ape::plot.phylo(x = data.tree,  # look at clustering
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
levels(gut) <- rev(c('12', '4', '10', '0', '19', '3', '7', '16', '18', '1', 
                     '11', '6', '2', '5', '9', '8', '13', '14', '15', '17')) 
gut$phyloorder <- Idents(gut) # create a new metadata slot with clusters in phylogenetic order

# Re-save Seurat object:
SaveH5Seurat(gut, filename = '/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpAll.h5Seurat', overwrite = TRUE)

# Plot some canonical genes:
Idents(gut) <- gut$phyloorder
DotPlot(gut, 
        features = c('EPCAM', 'CLDN3', 'KRT8', 'KRT20', 'CDH1',
                     'CD3E', 'CD3G', 'CD4', 'CD8B', 'CD8A', 'TRDC', 'CD2',
                     'CD79A', 'CD79B', 'CD19', 'MS4A1', 'JCHAIN', 
                     'CD14', 'SIRPA', 'MS4A2',
                     'PECAM1', 'COL3A1'),
        cols = c('yellow', 'darkgreen')) + RotatedAxis()
FeaturePlot(gut, 
            features = c('EPCAM', 'CLDN3', 'KRT8', 'KRT20', 'CDH1',
                         'CD3E', 'CD3G', 'CD4', 'CD8B', 'CD8A', 'TRDC', 'CD2',
                         'CD79A', 'CD79B', 'CD19', 'MS4A1', 'JCHAIN', 
                         'CD14', 'SIRPA', 'MS4A2',
                         'PECAM1', 'COL3A1'),
            cols = c('grey90', 'darkgreen'),
            ncol = 5,
            pt.size = 0.01) & NoLegend() & NoAxes()

# We see clusters 8, 9, 13, 14, 15, & 17 have gene expression characteristic of epithelial cells, but the occasional cell also expresses immune cell-specific genes, such as CD79B or TRDC
# We will identify epithelial cells as cells in clusters 8, 9, 13, 14, 15, & 17 that also have zero values for expression of selected leukocyte-specific genes

# Subset to only epithelial cells:
ep <- subset(gut, 
             idents = c('9', '8', '13', '14', '15', '17'), 
             subset = (CD3E == 0 & CD3G == 0 & CD4 == 0 & CD8B == 0 & CD8A == 0 & TRDC == 0 & CD2 == 0 & CD79A == 0 & CD79B == 0 & CD19 == 0 & MS4A1 == 0 & JCHAIN == 0 & SIRPA == 0 & CD14 == 0 & MS4A2 == 0  & PECAM1 == 0 & COL3A1 == 0))
ncol(ep) # how many cells left?
#[1] 695

# Point out only epithelial cells in our original/full dataset:
DimPlot(gut,
        cells.highlight = colnames(ep),
        cols.highlight = 'black')
DimPlot(gut,
        cells.highlight = colnames(ep),
        cols.highlight = 'black',
        split.by = 'orig.ident',
        ncol = 2) & NoLegend() & NoAxes()

# Now let's re-perform our normalization, integration, and dimensionality reduction on only our epithelial cells:
# Identify genes with non-zero expression to keep:
counts <- as.data.frame(ep[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])

# Slim down the Seurat object:
ep <- DietSeurat(ep, 
                  counts = TRUE,
                  data = TRUE,
                  scale.data = FALSE, # remove the scaled data
                  dimreducs = NULL,
                  features = keep, # keep only genes with non-zero counts across all cells
                  assays = 'RNA') # keep only RNA assay and remove SCT and integrated

# Re-normalize & re-integerate data:
ep.list <- SplitObject(ep, split.by = "orig.ident") # split by sample IDs
for (i in 1:length(ep.list)) { # normalize data using SCTransform method
  ep.list[[i]] <- SCTransform(ep.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}
ep.features <- SelectIntegrationFeatures(ep.list, # select the genes to use for integration
                                          verbose = TRUE) 
ep.list <- PrepSCTIntegration(ep.list, 
                               anchor.features = ep.features,
                               verbose = TRUE)
ep.anchors <- FindIntegrationAnchors(ep.list, # identify anchors for integration from top 30 PCs
                                      normalization.method = "SCT", 
                                      anchor.features = ep.features, 
                                      dims = 1:30,
                                     k.filter = 132) # reduce k.filter to number of cells in our smallest sample
ep.integrated <- IntegrateData(ep.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)

# Re-run PCA analysis:
ep.integrated <- RunPCA(ep.integrated, # run PCA analysis for 100 dimensions of the data
                         npcs = 100, 
                         verbose = TRUE) 
ElbowPlot(ep.integrated,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering
pct <- ep.integrated[["pca"]]@stdev / sum(ep.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
#co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
#co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
#[1] 10
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs 

# Re-do dimensionality reduction:
ep.integrated <- RunUMAP(ep.integrated, 
                          dims = PCdims, 
                          reduction = "pca", 
                          assay = "SCT") # create UMAP
ep.integrated <- RunTSNE(ep.integrated, 
                          dims = PCdims, 
                          reduction = "pca", 
                          assay = "SCT") # create tSNE plot (if desired)

# Re-cluster data:
ep.integrated <- FindNeighbors(ep.integrated, 
                                dims = PCdims, 
                                verbose = TRUE) 
ep.integrated <- FindClusters(ep.integrated, dims = PCdims, reduction = 'pca', resolution = 0.1) 

# Visualize clusters in UMAP & t-SNE
DimPlot(ep.integrated,
        group.by = 'seurat_clusters', 
        reduction = 'umap',
        label = TRUE)
DimPlot(ep.integrated,
        group.by = 'seurat_clusters', 
        reduction = 'tsne',
        label = TRUE)

# Visualize sample origins in UMAP & t-SNE
DimPlot(ep.integrated,
        group.by = 'orig.ident', 
        reduction = 'tsne')
DimPlot(ep.integrated,
        group.by = 'orig.ident', 
        reduction = 'tsne',
        cols = c('darkorange', 'dodgerblue3', 'chartreuse3', 'dodgerblue'))

# Re-normalize & re-scale RNA assay counts:
ep.integrated <- NormalizeData(ep.integrated,  
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
ep.integrated <- ScaleData(ep.integrated, 
                            assay = "RNA")

# Perform hierarchical clustering of cell clusters:
DefaultAssay(ep.integrated) <- 'RNA'
ep.integrated <- BuildClusterTree(ep.integrated, 
                        dims = PCdims, 
                        assay = "PCA")
data.tree <- Tool(object =ep.integrated, 
                  slot = "BuildClusterTree") 
ape::plot.phylo(x = data.tree,  # look at clustering
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)
data.tree <- ape::rotateConstr(data.tree, c('3', '0', '5', '1', '2', '4'))
plot(data.tree, direction = 'downwards', edge.width = 1.5, font = 1)
levels(ep.integrated) <- c('3', '0', '5', '1', '2', '4')
ep.integrated$phyloorder <- Idents(ep.integrated) # create a new metadata slot with clusters in phylogenetic order

SaveH5Seurat(ep.integrated, filename = "/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat")

#sessionInfo()

#R version 4.1.1 (2021-08-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggplot2_3.3.5         SeuratDisk_0.0.0.9019 SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#[1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              splines_4.1.1               BiocParallel_1.26.2         listenv_0.8.0              
#[7] scattermore_0.7             GenomeInfoDb_1.28.4         digest_0.6.27               htmltools_0.5.2             viridis_0.6.1               fansi_0.5.0                
#[13] magrittr_2.0.1              ScaledMatrix_1.0.0          tensor_1.5                  cluster_2.1.2               ROCR_1.0-11                 limma_3.48.3               
#[19] globals_0.14.0              graphlayouts_0.7.1          matrixStats_0.60.1          spatstat.sparse_2.0-0       colorspace_2.0-2            ggrepel_0.9.1              
#[25] xfun_0.26                   dplyr_1.0.7                 crayon_1.4.1                RCurl_1.98-1.4              jsonlite_1.7.2              spatstat.data_2.1-0        
#[31] ape_5.5                     survival_3.2-13             zoo_1.8-9                   glue_1.4.2                  polyclip_1.10-0             gtable_0.3.0               
#[37] zlibbioc_1.38.0             XVector_0.32.0              leiden_0.3.9                DelayedArray_0.18.0         BiocSingular_1.8.1          future.apply_1.8.1         
#[43] SingleCellExperiment_1.14.1 BiocGenerics_0.38.0         abind_1.4-5                 scales_1.1.1                edgeR_3.34.1                DBI_1.1.1                  
#[49] miniUI_0.1.1.1              Rcpp_1.0.7                  viridisLite_0.4.0           xtable_1.8-4                reticulate_1.21             spatstat.core_2.3-0        
#[55] rsvd_1.0.5                  bit_4.0.4                   htmlwidgets_1.5.4           httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.2             
#[61] ica_1.0-2                   pkgconfig_2.0.3             farver_2.1.0                uwot_0.1.10                 dbplyr_2.1.1                deldir_0.2-10              
#[67] locfit_1.5-9.4              utf8_1.2.2                  labeling_0.4.2              tidyselect_1.1.1            rlang_0.4.11                reshape2_1.4.4             
#[73] later_1.3.0                 munsell_0.5.0               cellranger_1.1.0            tools_4.1.1                 cli_3.0.1                   generics_0.1.0             
#[79] ggridges_0.5.3              evaluate_0.14               stringr_1.4.0               fastmap_1.1.0               yaml_2.2.1                  goftest_1.2-2              
#[85] knitr_1.34                  bit64_4.0.5                 fitdistrplus_1.1-5          tidygraph_1.2.0             purrr_0.3.4                 RANN_2.6.1                 
#[91] ggraph_2.0.5                pbapply_1.5-0               future_1.22.1               nlme_3.1-152                mime_0.11                   rstudioapi_0.13            
#[97] hdf5r_1.3.4                 compiler_4.1.1              beeswarm_0.4.0              plotly_4.9.4.1              png_0.1-7                   spatstat.utils_2.2-0       
#[103] tibble_3.1.4                tweenr_1.0.2                stringi_1.7.4               RSpectra_0.16-0             lattice_0.20-44             Matrix_1.3-4               
#[109] vctrs_0.3.8                 pillar_1.6.2                lifecycle_1.0.0             spatstat.geom_2.2-2         lmtest_0.9-38               BiocNeighbors_1.10.0       
#[115] RcppAnnoy_0.0.19            data.table_1.14.0           cowplot_1.1.1               bitops_1.0-7                irlba_2.3.3                 httpuv_1.6.3               
#[121] patchwork_1.1.1             GenomicRanges_1.44.0        R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3              
#[127] vipor_0.4.5                 IRanges_2.26.0              parallelly_1.28.1           codetools_0.2-18            gtools_3.9.2                MASS_7.3-54                
#[133] assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.4.2                 sctransform_0.3.2           S4Vectors_0.30.0            GenomeInfoDbData_1.2.6     
#[139] mgcv_1.8-36                 miloR_1.0.0                 beachmat_2.8.1              grid_4.1.1                  rpart_4.1-15                tidyr_1.1.3                
#[145] rmarkdown_2.11              MatrixGenerics_1.4.3        Rtsne_0.15                  ggforce_0.3.3               Biobase_2.52.0              shiny_1.6.0                
#[151] ggbeeswarm_0.6.0     
