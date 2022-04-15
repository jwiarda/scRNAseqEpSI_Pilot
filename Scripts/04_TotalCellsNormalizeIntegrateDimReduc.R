# Title: 04 - Normalization, Integration, Dimensionality Reduction of Total Cells
# Date: 2021Sept21
# Author: Jayne Wiarda

library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)

## Import & split Seurat object:
All <- readRDS("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/GutEpScrubbedSeurat.rds") # read in Seurat object from doublet removal
All.list <- SplitObject(All, split.by = "orig.ident") # split by sample IDs
All.list <- All.list[c("DUOD2", "JEJ2", "IPP2", "NoPP2")] 

## Perform SCTransform normalization of data from each sample:
for (i in 1:length(All.list)) { # normalize data using SCTransform method
  All.list[[i]] <- SCTransform(All.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}

## Integrate the data from different samples:
All.features <- SelectIntegrationFeatures(All.list, # select the genes to use for integration
                                          verbose = TRUE) 
All.list <- PrepSCTIntegration(All.list, 
                               anchor.features = All.features,
                               verbose = TRUE)
All.anchors <- FindIntegrationAnchors(All.list, # identify anchors for integration from top 30 PCs
                                      normalization.method = "SCT", 
                                      anchor.features = All.features, 
                                      dims = 1:30)
All.integrated <- IntegrateData(All.anchors, # integrate data
                                normalization.method = "SCT", 
                                dims = 1:30)

## Run multidimensional analyses on data:
All.integrated <- RunPCA(All.integrated, # run PCA analysis for 50 dimensions of the data
                         npcs = 100, 
                         verbose = TRUE) 
ElbowPlot(All.integrated,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering

## Quantitiatively calculate your PC:
pct <- All.integrated[["pca"]]@stdev / sum(All.integrated[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co1 # list PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
co2 # list PC
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
pcs # list PC
plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps

## Perform multidimensional visualization of data:
All.integrated <- RunUMAP(All.integrated, 
                          dims = PCdims, 
                          reduction = "pca", 
                          assay = "SCT") # create UMAP
All.integrated <- RunTSNE(All.integrated, 
                          dims = PCdims, 
                          reduction = "pca", 
                          assay = "SCT") # create tSNE plot 

## Also define nearest neighbors:
All.integrated <- FindNeighbors(All.integrated, 
                                dims = PCdims, 
                                verbose = TRUE) 

## Add normalized/scaled data to RNA assay:
#dim(All.integrated[["RNA"]]@scale.data) # see that there is no RNA assay scaled data yet
All.integrated <- NormalizeData(All.integrated,  # normalize the RNA counts data per cell
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000, 
                                assay = "RNA")
All.integrated <- ScaleData(All.integrated, # scale the RNA counts data relative to other cells
                            assay = "RNA")
#dim(All.integrated[["RNA"]]@scale.data) # see that all genes are scaled in RNA assay now

DefaultAssay(All.integrated) <- "RNA"

# Visualize dimensionality reductions of integrated data:

DimPlot(All.integrated, label = FALSE, group.by = "orig.ident")
DimPlot(All.integrated, label = FALSE, split.by = "orig.ident")
DimPlot(All.integrated, label = FALSE, reduction = 'tsne', group.by = "orig.ident")
DimPlot(All.integrated, label = FALSE, reduction = 'tsne', split.by = "orig.ident")

## Save the Seurat object:
SaveH5Seurat(All.integrated, filename = "/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpAll.h5Seurat")

#sessionInfo()
#R version 4.1.1 (2021-08-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggplot2_3.3.5         writexl_1.4.0         dplyr_1.0.7           SeuratDisk_0.0.0.9019 SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#[1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              splines_4.1.1               BiocParallel_1.26.2         listenv_0.8.0               scattermore_0.7             GenomeInfoDb_1.28.4         digest_0.6.27              
#[10] htmltools_0.5.2             viridis_0.6.1               fansi_0.5.0                 magrittr_2.0.1              ScaledMatrix_1.0.0          tensor_1.5                  cluster_2.1.2               ROCR_1.0-11                 limma_3.48.3               
#[19] globals_0.14.0              graphlayouts_0.7.1          matrixStats_0.60.1          spatstat.sparse_2.0-0       colorspace_2.0-2            ggrepel_0.9.1               xfun_0.26                   crayon_1.4.1                RCurl_1.98-1.4             
#[28] jsonlite_1.7.2              spatstat.data_2.1-0         survival_3.2-13             zoo_1.8-9                   glue_1.4.2                  polyclip_1.10-0             gtable_0.3.0                zlibbioc_1.38.0             XVector_0.32.0             
#[37] leiden_0.3.9                DelayedArray_0.18.0         BiocSingular_1.8.1          future.apply_1.8.1          SingleCellExperiment_1.14.1 BiocGenerics_0.38.0         abind_1.4-5                 scales_1.1.1                edgeR_3.34.1               
#[46] DBI_1.1.1                   miniUI_0.1.1.1              Rcpp_1.0.7                  viridisLite_0.4.0           xtable_1.8-4                reticulate_1.21             spatstat.core_2.3-0         rsvd_1.0.5                  bit_4.0.4                  
#[55] htmlwidgets_1.5.4           httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.2              ica_1.0-2                   pkgconfig_2.0.3             farver_2.1.0                uwot_0.1.10                 dbplyr_2.1.1               
#[64] deldir_0.2-10               locfit_1.5-9.4              utf8_1.2.2                  labeling_0.4.2              tidyselect_1.1.1            rlang_0.4.11                reshape2_1.4.4              later_1.3.0                 munsell_0.5.0              
#[73] cellranger_1.1.0            tools_4.1.1                 cli_3.0.1                   generics_0.1.0              ggridges_0.5.3              evaluate_0.14               stringr_1.4.0               fastmap_1.1.0               yaml_2.2.1                 
#[82] goftest_1.2-2               knitr_1.34                  bit64_4.0.5                 fitdistrplus_1.1-5          tidygraph_1.2.0             purrr_0.3.4                 RANN_2.6.1                  ggraph_2.0.5                pbapply_1.5-0              
#[91] future_1.22.1               nlme_3.1-152                mime_0.11                   hdf5r_1.3.4                 compiler_4.1.1              beeswarm_0.4.0              plotly_4.9.4.1              png_0.1-7                   spatstat.utils_2.2-0       
#[100] tibble_3.1.4                tweenr_1.0.2                stringi_1.7.4               lattice_0.20-44             Matrix_1.3-4                vctrs_0.3.8                 pillar_1.6.2                lifecycle_1.0.0             spatstat.geom_2.2-2        
#[109] lmtest_0.9-38               RcppAnnoy_0.0.19            BiocNeighbors_1.10.0        data.table_1.14.0           cowplot_1.1.1               bitops_1.0-7                irlba_2.3.3                 httpuv_1.6.3                patchwork_1.1.1            
#[118] GenomicRanges_1.44.0        R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3               vipor_0.4.5                 IRanges_2.26.0              parallelly_1.28.1           codetools_0.2-18           
#[127] gtools_3.9.2                MASS_7.3-54                 assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.4.2                 sctransform_0.3.2           S4Vectors_0.30.0            GenomeInfoDbData_1.2.6      mgcv_1.8-36                
#[136] miloR_1.0.0                 beachmat_2.8.1              grid_4.1.1                  rpart_4.1-15                tidyr_1.1.3                 rmarkdown_2.11              MatrixGenerics_1.4.3        Rtsne_0.15                  ggforce_0.3.3              
#[145] Biobase_2.52.0              shiny_1.6.0                 ggbeeswarm_0.6.0      