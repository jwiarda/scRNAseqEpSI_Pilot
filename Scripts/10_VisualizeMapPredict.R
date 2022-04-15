# Title: 10 - Visualizing Epithelial Cell Mapping & Prediction Results
# Date: 2021Nov23
# Author: Jayne Wiarda

library(Seurat)
library(SeuratDisk)
library(readxl)
library(ggplot2)
library(scales)
library(tidyr)
library(viridis)

# Load in dataset:
ep <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
DefaultAssay <- 'RNA'
Idents(ep) <- ep$celltype

# Incorporate mouse mapping/prediction information:
MappingScores <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/MmSI_EpOnly_Haber2017_MappingScores.xlsx')
CellTypePredictions <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/MmSI_EpOnly_Haber2017_CellTypePredictions.xlsx')
CellTypePredictions$prediction.score.enterocyte.sum <- rowSums(CellTypePredictions[ , c('prediction.score.EP', 'prediction.score.Enterocyte')], na.rm=TRUE) # sum all enterocyte subtypes into a single enterocyte prediction
CellTypePredictions$prediction.score.crypt.sum <- rowSums(CellTypePredictions[ , c('prediction.score.Stem', 'prediction.score.TA', 'prediction.score.Paneth')], na.rm=TRUE) # sum all crypt cell subtypes into a single crypt cell prediction
colnames(MappingScores) <-paste("Mm", colnames(MappingScores), sep="_")
colnames(CellTypePredictions) <-paste("Mm", colnames(CellTypePredictions), sep="_")

# Incorporate cell prediction & mapping scores into original Seurat object of query data:
ep <- AddMetaData(object = ep, 
                  metadata = c(MappingScores, CellTypePredictions))

# Go back through and re-calculate predicted id based on our merged prediction score amendments:
colnames(CellTypePredictions)
CellBarcodes <- CellTypePredictions$Mm_CellBarcodes
CellTypePredictions <- subset(CellTypePredictions, select=-c(Mm_predicted.id, Mm_prediction.score.max, Mm_CellBarcodes,
                                                             Mm_prediction.score.EP, Mm_prediction.score.Enterocyte,
                                                             Mm_prediction.score.Stem, Mm_prediction.score.TA, Mm_prediction.score.Paneth))
colnames(CellTypePredictions)
RevisedPredictedIDs <- as.data.frame(colnames(CellTypePredictions[max.col(CellTypePredictions, ties.method = "first")]))
RevisedPredictedIDs$CellBarcodes <- CellBarcodes
colnames(RevisedPredictedIDs) <- c('MmPredictedIdRevised', 'CellBarcodes')
ep <- AddMetaData(object = ep, 
                  metadata = c(RevisedPredictedIDs))

# Incorporate human mapping/prediction information:
MappingScores <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/HsSI_EpOnly_Elmentaite2021_MappingScores.xlsx')
CellTypePredictions <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/HsSI_EpOnly_Elmentaite2021_CellTypePredictions.xlsx')
CellTypePredictions$prediction.score.EE.sum <- rowSums(CellTypePredictions[ , c('prediction.score.L.cells..PYY..', 'prediction.score.EC.cells..TAC1..',
                                                                                'prediction.score.D.cells..SST..', 'prediction.score.N.cells..NTS..',
                                                                                'prediction.score.Progenitor..NEUROG3..', 'prediction.score.I.cells..CCK..',
                                                                                'prediction.score.K.cells..GIP..', 'prediction.score.M.X.cells..MLN.GHRL..',
                                                                                'prediction.score.EECs')], na.rm=TRUE) # sum all EE cell subtypes into a single EE cell prediction
CellTypePredictions$prediction.score.goblet.sum <- rowSums(CellTypePredictions[ , c('prediction.score.BEST2..Goblet.cell', 'prediction.score.Goblet.cell')], na.rm=TRUE) # sum all goblet cell subtypes into a single goblet cell prediction
CellTypePredictions$prediction.score.crypt.sum <- rowSums(CellTypePredictions[ , c('prediction.score.TA', 'prediction.score.Stem.cells', 
                                                                                   'prediction.score.Paneth')], na.rm=TRUE) # sum all crypt cell subtypes into a single crypt cell prediction
colnames(MappingScores) <-paste("Hs", colnames(MappingScores), sep="_")
colnames(CellTypePredictions) <-paste("Hs", colnames(CellTypePredictions), sep="_")

# Incorporate cell prediction & mapping scores into original Seurat object of query data:
ep <- AddMetaData(object = ep, 
                  metadata = c(MappingScores, CellTypePredictions))

# Go back through and re-calculate predicted id based on our merged prediction score amendments:
colnames(CellTypePredictions)
CellBarcodes <- CellTypePredictions$Hs_CellBarcodes
CellTypePredictions <- subset(CellTypePredictions, select=-c(Hs_predicted.id, Hs_prediction.score.L.cells..PYY.., Hs_prediction.score.EC.cells..TAC1..,
                                                             Hs_prediction.score.D.cells..SST.., Hs_prediction.score.N.cells..NTS..,
                                                             Hs_prediction.score.Progenitor..NEUROG3.., Hs_prediction.score.I.cells..CCK..,
                                                             Hs_prediction.score.K.cells..GIP.., Hs_prediction.score.M.X.cells..MLN.GHRL..,
                                                             Hs_prediction.score.EECs,
                                                             Hs_prediction.score.BEST2..Goblet.cell, Hs_prediction.score.Goblet.cell,
                                                             Hs_prediction.score.TA, Hs_prediction.score.Stem.cells, 
                                                             Hs_prediction.score.Paneth, Hs_prediction.score.max, Hs_CellBarcodes))
colnames(CellTypePredictions)
RevisedPredictedIDs <- as.data.frame(colnames(CellTypePredictions[max.col(CellTypePredictions, ties.method = "first")]))
RevisedPredictedIDs$CellBarcodes <- CellBarcodes
colnames(RevisedPredictedIDs) <- c('HsPredictedIdRevised', 'CellBarcodes')
ep <- AddMetaData(object = ep, 
                  metadata = c(RevisedPredictedIDs))

# Plot mapping scores:
FeaturePlot(ep, features = c('Hs_MappingScores', 'Mm_MappingScores'), 
            reduction = 'tsne', 
            cols = c('khaki1', 'yellow', 'gold', 'orange', 'red', 'darkred'),
            min.cutoff = 0, 
            max.cutoff = 1) & NoAxes() & DarkTheme()
VlnPlot(ep, features = c('Hs_MappingScores', 'Mm_MappingScores'), 
        cols = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta'),
        y.max = 1,
        pt.size = 0.01) & 
  geom_boxplot(width = 0.15, fill = 'white', col = 'black', lwd = 0.5, outlier.shape = NA, coef = 0) &
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "red") # red dot shows mean; box shows IQR

# Plot predicted cell ID (based on highest prediction scores):
DimPlot(ep, 
        group.by = 'HsPredictedIdRevised',
        reduction = 'tsne',
        cols = c('darkorange2', 'burlywood4', 'blue', 'goldenrod1', 'chartreuse3', 'red', 'grey50'))
DimPlot(ep, 
        group.by = 'MmPredictedIdRevised',
        reduction = 'tsne',
        cols = c('burlywood4', 'goldenrod1', 'blue', 'chartreuse3', 'grey50'))

# Plot prediction scores for different cell types from reference datasets:
FeaturePlot(ep,
            features = "Mm_prediction.score.crypt.sum",
            reduction = 'tsne') & 
  scale_color_gradientn( colours = c('grey90', 'darkslateblue'),  
                         limits = c(0, 1), 
                         oob = squish) & 
  NoAxes() & NoLegend() 
VlnPlot(ep, features = c('Mm_prediction.score.crypt.sum'), 
        cols = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta'),
        y.max = 1,
        pt.size = 0.01) & 
  geom_boxplot(width = 0.15, fill = 'white', col = 'black', lwd = 0.5, outlier.shape = NA, coef = 0) &
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "red")

# Previous two commands can be repeated with substitution of other cell types for plotting. Other available cell types listed below:
"Hs_prediction.score.Enterocyte"
"Hs_prediction.score.TA"                    
"Hs_prediction.score.Stem.cells"           
"Hs_prediction.score.Goblet.cell"           
"Hs_prediction.score.L.cells..PYY.."        
"Hs_prediction.score.Tuft"                 
"Hs_prediction.score.Microfold.cell"        
"Hs_prediction.score.EC.cells..TAC1.."      
"Hs_prediction.score.D.cells..SST.."       
"Hs_prediction.score.N.cells..NTS.."        
"Hs_prediction.score.BEST4..epithelial"     
"Hs_prediction.score.Paneth"               
"Hs_prediction.score.Progenitor..NEUROG3.." 
"Hs_prediction.score.BEST2..Goblet.cell"    
"Hs_prediction.score.EECs"                 
"Hs_prediction.score.I.cells..CCK.."      
"Hs_prediction.score.K.cells..GIP.."        
"Hs_prediction.score.M.X.cells..MLN.GHRL.."
"Hs_prediction.score.EE.sum"                
"Hs_prediction.score.crypt.sum"            
"Hs_prediction.score.goblet.sum"    
"Mm_prediction.score.Goblet"          
"Mm_prediction.score.Enterocyte"      
"Mm_prediction.score.Stem"           
"Mm_prediction.score.TA"              
"Mm_prediction.score.Paneth"          
"Mm_prediction.score.EP"             
"Mm_prediction.score.Tuft"            
"Mm_prediction.score.Enteroendocrine"         
"Mm_prediction.score.enterocyte.sum"  
"Mm_prediction.score.crypt.sum"   

# Create heat map of average prediction scores:
meta <- ep@meta.data
metaMm <- subset(meta, select = c("celltype", "Mm_prediction.score.Goblet", "Mm_prediction.score.Tuft", "Mm_prediction.score.Enteroendocrine",
                                  "Mm_prediction.score.enterocyte.sum", "Mm_prediction.score.crypt.sum"))
metaMm <- aggregate(metaMm, by = list(metaMm$celltype), mean)
metaMm <- gather(metaMm, key = "refCellType", value = "meanPrediction",
                 Mm_prediction.score.Goblet, Mm_prediction.score.Tuft, Mm_prediction.score.Enteroendocrine,
                                  Mm_prediction.score.enterocyte.sum, Mm_prediction.score.crypt.sum)
metaMm$refCellType <- factor(metaMm$refCellType, levels = c("Mm_prediction.score.crypt.sum", "Mm_prediction.score.enterocyte.sum", 
                                                            "Mm_prediction.score.Goblet", "Mm_prediction.score.Enteroendocrine", 
                                                            "Mm_prediction.score.Tuft"))
ggplot(metaMm, aes(refCellType, Group.1, fill= meanPrediction)) + 
  geom_tile(color = 'grey30', lwd = 0.5) +
  scale_fill_viridis(option = 'viridis', limits = c(0, 0.7), oob = squish) + 
  theme_classic() + 
  RotatedAxis()

metaHs <- subset(meta, select = c("celltype", "Hs_prediction.score.Enterocyte", "Hs_prediction.score.Tuft",              
                                  "Hs_prediction.score.Microfold.cell", "Hs_prediction.score.BEST4..epithelial",
                                  "Hs_prediction.score.EE.sum", "Hs_prediction.score.goblet.sum", "Hs_prediction.score.crypt.sum"))
metaHs <- aggregate(metaHs, by = list(metaHs$celltype), mean)
metaHs <- gather(metaHs, key = "refCellType", value = "meanPrediction",
                 Hs_prediction.score.Enterocyte, Hs_prediction.score.Tuft,              
                 Hs_prediction.score.Microfold.cell, Hs_prediction.score.BEST4..epithelial,
                 Hs_prediction.score.EE.sum, Hs_prediction.score.goblet.sum, Hs_prediction.score.crypt.sum)
metaHs$refCellType <- factor(metaHs$refCellType, levels = c("Hs_prediction.score.crypt.sum", "Hs_prediction.score.Enterocyte", 
                                                    "Hs_prediction.score.BEST4..epithelial", "Hs_prediction.score.goblet.sum",
                                                    "Hs_prediction.score.EE.sum", "Hs_prediction.score.Tuft", "Hs_prediction.score.Microfold.cell"))
ggplot(metaHs, aes(refCellType, Group.1, fill= meanPrediction)) + 
  geom_tile(color = 'grey30', lwd = 0.5) +
  scale_fill_viridis(option = 'viridis', limits = c(0, 0.7), oob = squish) + 
  theme_classic() + 
  RotatedAxis()

sessionInfo()
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] viridis_0.6.2         viridisLite_0.4.0     tidyr_1.1.3           scales_1.1.1          ggplot2_3.3.5         readxl_1.3.1          SeuratDisk_0.0.0.9019
#[8] SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#[1] scattermore_0.7             pkgmaker_0.32.2             bit64_4.0.5                 knitr_1.34                  multcomp_1.4-17             irlba_2.3.3                
#[7] DelayedArray_0.18.0         data.table_1.14.0           rpart_4.1-15                KEGGREST_1.32.0             RCurl_1.98-1.4              doParallel_1.0.16          
#[13] generics_0.1.0              metap_1.5                   BiocGenerics_0.38.0         ScaledMatrix_1.0.0          TH.data_1.1-0               cowplot_1.1.1              
#[19] RSQLite_2.2.8               RANN_2.6.1                  miloR_1.0.0                 future_1.22.1               bit_4.0.4                   mutoss_0.1-12              
#[25] phylobase_0.8.10            spatstat.data_2.1-0         xml2_1.3.2                  httpuv_1.6.3                SummarizedExperiment_1.22.0 assertthat_0.2.1           
#[31] xfun_0.26                   hms_1.1.0                   evaluate_0.14               promises_1.2.0.1            fansi_0.5.0                 progress_1.2.2             
#[37] igraph_1.2.6                DBI_1.1.1                   tmvnsim_1.0-2               htmlwidgets_1.5.4           spatstat.geom_2.2-2         stats4_4.1.2               
#[43] purrr_0.3.4                 ellipsis_0.3.2              dplyr_1.0.7                 annotate_1.70.0             gridBase_0.4-7              locfdr_1.1-8               
#[49] deldir_0.2-10               MatrixGenerics_1.4.3        vctrs_0.3.8                 SingleCellExperiment_1.14.1 Biobase_2.52.0              ROCR_1.0-11                
#[55] abind_1.4-5                 cachem_1.0.6                withr_2.4.2                 ggforce_0.3.3               sctransform_0.3.2           prettyunits_1.1.1          
#[61] goftest_1.2-2               mnormt_2.0.2                softImpute_1.4-1            cluster_2.1.2               ape_5.5                     lazyeval_0.2.2             
#[67] crayon_1.4.1                genefilter_1.74.0           hdf5r_1.3.4                 edgeR_3.34.1                pkgconfig_2.0.3             labeling_0.4.2             
#[73] tweenr_1.0.2                GenomeInfoDb_1.28.4         nlme_3.1-152                vipor_0.4.5                 rlang_0.4.11                globals_0.14.0             
#[79] lifecycle_1.0.0             miniUI_0.1.1.1              sandwich_3.0-1              registry_0.5-1              mathjaxr_1.4-0              rsvd_1.0.5                 
#[85] cellranger_1.1.0            polyclip_1.10-0             matrixStats_0.60.1          lmtest_0.9-38               rngtools_1.5                Matrix_1.3-4               
#[91] Rhdf5lib_1.14.2             zoo_1.8-9                   beeswarm_0.4.0              ggridges_0.5.3              png_0.1-7                   bitops_1.0-7               
#[97] rncl_0.8.4                  KernSmooth_2.23-20          rhdf5filters_1.4.0          Biostrings_2.60.2           blob_1.2.2                  stringr_1.4.0              
#[103] zinbwave_1.14.2             parallelly_1.28.1           S4Vectors_0.30.2            beachmat_2.8.1              memoise_2.0.0               magrittr_2.0.1             
#[109] plyr_1.8.6                  ica_1.0-2                   howmany_0.3-1               zlibbioc_1.38.0             compiler_4.1.2              RColorBrewer_1.1-2         
#[115] plotrix_3.8-2               fitdistrplus_1.1-5          cli_3.0.1                   ade4_1.7-18                 XVector_0.32.0              listenv_0.8.0              
#[121] patchwork_1.1.1             pbapply_1.5-0               MASS_7.3-54                 mgcv_1.8-38                 tidyselect_1.1.1            stringi_1.7.4              
#[127] yaml_2.2.1                  BiocSingular_1.8.1          locfit_1.5-9.4              ggrepel_0.9.1               grid_4.1.2                  tools_4.1.2                
#[133] future.apply_1.8.1          rstudioapi_0.13             uuid_0.1-4                  foreach_1.5.1               RNeXML_2.4.5                gridExtra_2.3              
#[139] farver_2.1.0                Rtsne_0.15                  ggraph_2.0.5                digest_0.6.27               shiny_1.6.0                 Rcpp_1.0.7                 
#[145] GenomicRanges_1.44.0        later_1.3.0                 RcppAnnoy_0.0.19            httr_1.4.2                  AnnotationDbi_1.54.1        kernlab_0.9-29             
#[151] Rdpack_2.1.2                colorspace_2.0-2            XML_3.99-0.7                tensor_1.5                  reticulate_1.21             clusterExperiment_2.12.0   
#[157] IRanges_2.26.0              splines_4.1.2               uwot_0.1.10                 sn_2.0.0                    spatstat.utils_2.2-0        graphlayouts_0.7.1         
#[163] multtest_2.48.0             plotly_4.9.4.1              xtable_1.8-4                jsonlite_1.7.2              tidygraph_1.2.0             R6_2.5.1                   
#[169] TFisher_0.2.0               pillar_1.6.2                htmltools_0.5.2             mime_0.11                   NMF_0.23.0                  glue_1.4.2                 
#[175] fastmap_1.1.0               BiocParallel_1.26.2         BiocNeighbors_1.10.0        codetools_0.2-18            mvtnorm_1.1-2               utf8_1.2.2                 
#[181] lattice_0.20-45             spatstat.sparse_2.0-0       tibble_3.1.4                numDeriv_2016.8-1.1         ggbeeswarm_0.6.0            leiden_0.3.9               
#[187] gtools_3.9.2                survival_3.2-13             limma_3.48.3                rmarkdown_2.11              munsell_0.5.0               rhdf5_2.36.0               
#[193] GenomeInfoDbData_1.2.6      iterators_1.0.13            HDF5Array_1.20.0            reshape2_1.4.4              gtable_0.3.0                rbibutils_2.2.3            
#[199] spatstat.core_2.3-0      