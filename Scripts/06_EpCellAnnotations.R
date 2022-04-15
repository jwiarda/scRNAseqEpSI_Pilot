# Title: 06 - Annotating Epithelial Cells
# Date: 2021Oct06
# Author: Jayne Wiarda

library(Seurat)  
library(SeuratDisk)

# Load in dataset:
ep <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
DefaultAssay <- 'RNA'
Idents(ep) <- ep$phyloorder

# After some gene query, we came up with annotations for our epithelial cell clusters:
ep$celltype <-ep$phyloorder
Idents(ep) <- ep$celltype
types <- c('Crypt', 'Enterocyte', 'BEST4 enterocyte', 'Goblet', 'NEUROD1lo EE', 'NEUROD1hi EE')
names(types) <- levels(ep) # assign GutCellTypes to cluster numbers
ep <- RenameIdents(ep, types) # change dataset identity to cell types in new Seurat object
ep$celltype <- Idents(ep)
Idents(ep) <- ep$celltype
DimPlot(ep,
        reduction = 'tsne',
        cols = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta'))
DimPlot(ep,
        reduction = 'tsne',
        cols = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta'),
        split.by = 'orig.ident') & NoAxes() & NoLegend()

# Also create meta data annotation for intestinal region:
ep$region <-ep$orig.ident
Idents(ep) <- ep$region
regions <- c('Duodenum', 'Jejunum', 'Ileum', 'Ileum')
names(regions) <- levels(ep) # assign GutCellTypes to cluster numbers
ep <- RenameIdents(ep, regions) # change dataset identity to cell types in new Seurat object
ep$region <- Idents(ep)
Idents(ep) <- ep$region
levels(ep) <- c('Duodenum', 'Jejunum', 'Ileum') 
ep$region <- Idents(ep)
Idents(ep) <- ep$region

# Split cells in UMAP by intestinal region:
Idents(ep) <- ep$celltype
DimPlot(ep,
        reduction = 'tsne',
        cols = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta'),
        split.by = 'region') & NoAxes() & NoLegend()

# Stacked barplots of cells derived from each intestinal region within each cell type:
Idents(ep) <- ep$celltype
GutregionTotalCells <- prop.table(table(ep$region)) # What percent of total cells are from each region?
GutregionPercents <- prop.table(table(Idents(ep),ep$region), 
                                margin = 1) # What percent of cells from each cluster belong to each region?
GutregionPercents <- rbind(GutregionPercents, GutregionTotalCells) # add row of overall percentages to table
#rowSums(GutregionPercents) # make sure all row sums are equal to 1
GutregionPercents <- t(GutregionPercents) # transpose the table
par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(GutregionPercents, # create stacked bar plot
        col = c('grey20', 'grey50', 'grey80'), 
        legend = rownames(GutregionPercents),
        ylab = "Frequency within cell type", 
        las = 2,
        border = NA,
        space = 0.05,
        legend.text = TRUE, 
        args.legend = list(x = "topright", bty = "n", inset=c(-0.5, 0)))

# Stacked bar of cell types within each original sample:
Idents(ep) <- ep$orig.ident
GutregionTotalCells <- prop.table(table(ep$celltype)) # What percent of total cells are from each region?
GutregionPercents <- prop.table(table(Idents(ep),ep$celltype), 
                                margin = 1) # What percent of cells from each cluster belong to each region?
GutregionPercents <- rbind(GutregionPercents, GutregionTotalCells) # add row of overall percentages to table
#rowSums(GutregionPercents) # make sure all row sums are equal to 1
GutregionPercents <- t(GutregionPercents) # transpose the table
par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(GutregionPercents, # create stacked bar plot
        col = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta'), 
        legend = rownames(GutregionPercents),
        ylab = "Frequency within cell type", 
        las = 2,
        border = NA,
        space = 0.05,
        legend.text = TRUE, 
        args.legend = list(x = "topright", bty = "n", inset=c(-0.6, 0)))

# Save over the h5Seurat object to incorporate new annotations:
SaveH5Seurat(ep, filename = "/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat", overwrite = TRUE) # overwrite = TRUE writes over previously-saved object, so use w/ caution

sessionInfo()
#R version 4.1.1 (2021-08-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] SeuratDisk_0.0.0.9019 SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#[1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              splines_4.1.1               BiocParallel_1.26.2        
#[6] listenv_0.8.0               scattermore_0.7             GenomeInfoDb_1.28.4         ggplot2_3.3.5               digest_0.6.27              
#[11] htmltools_0.5.2             viridis_0.6.1               fansi_0.5.0                 magrittr_2.0.1              ScaledMatrix_1.0.0         
#[16] tensor_1.5                  cluster_2.1.2               ROCR_1.0-11                 limma_3.48.3                globals_0.14.0             
#[21] graphlayouts_0.7.1          matrixStats_0.60.1          spatstat.sparse_2.0-0       colorspace_2.0-2            ggrepel_0.9.1              
#[26] xfun_0.26                   dplyr_1.0.7                 crayon_1.4.1                RCurl_1.98-1.4              jsonlite_1.7.2             
#[31] spatstat.data_2.1-0         ape_5.5                     survival_3.2-13             zoo_1.8-9                   glue_1.4.2                 
#[36] polyclip_1.10-0             gtable_0.3.0                zlibbioc_1.38.0             XVector_0.32.0              leiden_0.3.9               
#[41] DelayedArray_0.18.0         BiocSingular_1.8.1          future.apply_1.8.1          SingleCellExperiment_1.14.1 BiocGenerics_0.38.0        
#[46] abind_1.4-5                 scales_1.1.1                edgeR_3.34.1                DBI_1.1.1                   miniUI_0.1.1.1             
#[51] Rcpp_1.0.7                  viridisLite_0.4.0           xtable_1.8-4                reticulate_1.21             spatstat.core_2.3-0        
#[56] rsvd_1.0.5                  bit_4.0.4                   htmlwidgets_1.5.4           httr_1.4.2                  RColorBrewer_1.1-2         
#[61] ellipsis_0.3.2              ica_1.0-2                   pkgconfig_2.0.3             farver_2.1.0                uwot_0.1.10                
#[66] dbplyr_2.1.1                deldir_0.2-10               locfit_1.5-9.4              utf8_1.2.2                  labeling_0.4.2             
#[71] tidyselect_1.1.1            rlang_0.4.11                reshape2_1.4.4              later_1.3.0                 munsell_0.5.0              
#[76] tools_4.1.1                 cli_3.0.1                   generics_0.1.0              ggridges_0.5.3              evaluate_0.14              
#[81] stringr_1.4.0               fastmap_1.1.0               yaml_2.2.1                  goftest_1.2-2               knitr_1.34                 
#[86] bit64_4.0.5                 fitdistrplus_1.1-5          tidygraph_1.2.0             purrr_0.3.4                 RANN_2.6.1                 
#[91] ggraph_2.0.5                pbapply_1.5-0               future_1.22.1               nlme_3.1-152                mime_0.11                  
#[96] hdf5r_1.3.4                 compiler_4.1.1              beeswarm_0.4.0              plotly_4.9.4.1              png_0.1-7                  
#[101] spatstat.utils_2.2-0        tibble_3.1.4                tweenr_1.0.2                stringi_1.7.4               lattice_0.20-45            
#[106] Matrix_1.3-4                vctrs_0.3.8                 pillar_1.6.2                lifecycle_1.0.0             spatstat.geom_2.2-2        
#[111] lmtest_0.9-38               RcppAnnoy_0.0.19            BiocNeighbors_1.10.0        data.table_1.14.0           cowplot_1.1.1              
#[116] bitops_1.0-7                irlba_2.3.3                 httpuv_1.6.3                patchwork_1.1.1             GenomicRanges_1.44.0       
#[121] R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3               vipor_0.4.5                
#[126] IRanges_2.26.0              parallelly_1.28.1           codetools_0.2-18            gtools_3.9.2                MASS_7.3-54                
#[131] assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.4.2                 sctransform_0.3.2           S4Vectors_0.30.2           
#[136] GenomeInfoDbData_1.2.6      mgcv_1.8-37                 miloR_1.0.0                 beachmat_2.8.1              grid_4.1.1                 
#[141] rpart_4.1-15                tidyr_1.1.3                 rmarkdown_2.11              MatrixGenerics_1.4.3        Rtsne_0.15                 
#[146] ggforce_0.3.3               Biobase_2.52.0              shiny_1.6.0                 ggbeeswarm_0.6.0   