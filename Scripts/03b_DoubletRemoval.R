# Title: 03b - Doublet Removal
# Date: 2021Sept21
# Author: Jayne Wiarda

library(dplyr)        
library(DropletUtils) 
library(Seurat)       

## Read in our .rds file from QC script and re-filter our counts and pheno data based on the previous pass/fail criteria. This leaves us with the same lists of cells we put into Scrublet.
All <- readRDS("/home/Jayne.Wiarda/scRNAseqEpSI/QC/GutEpQC.rds")
#dim(All$counts)
All$counts <- All$counts[,All$phenoData$PassAll] # keep only cells that passed all criteria in QC
#dim(All$counts) 
All$phenoData <- All$phenoData[All$phenoData$PassAll,]
#dim(All$phenoData) # make sure rows here matches columns in All$counts

## Import our Scrublet scores:
scrubDUOD2 <- read.csv("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/DUOD2_ScrubScore.csv") # import doublet scores for sample DUOD2
scrubJEJ2 <- read.csv("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/JEJ2_ScrubScore.csv") # import doublet scores for sample JEJ2
scrubIPP2 <- read.csv("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/IPP2_ScrubScore.csv") # import doublet scores for sample IPP2
scrubNoPP2 <- read.csv("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/NoPP2_ScrubScore.csv") # import doublet scores for sample NOPP2
scrub_combined <- rbind(scrubDUOD2, scrubJEJ2, scrubIPP2, scrubNoPP2) # bind the dataframes together in the exact same order we've been using to list our samples in previous scripts
#dim(scrub_combined) # check the dimensions for number of rows is identical to the number of rows in All$phenoData and columns in All$counts
#head(scrub_combined, n = 3) # See what output looks like to identify column of ScrubScores. This column is called X0

All$phenoData$Scrublet <- scrub_combined$X0 # add column of scrublet data information for each cell to phenotype dataframe
All$phenoData <- mutate(All$phenoData, PassScrub = Scrublet < 0.25) # set the appropriate probability threshold... we decided on 0.25 for all samples here
rownames(All$phenoData) <- All$phenoData$Loupe # make sure the All$phenoData rownames correspond to Loupe barcodes!

#table(All$phenoData$SampleID,All$phenoData$PassScrub) # how many cells passed scrubbing criteria?

## Filter out probable doublets:
#dim(All$counts)
All$counts <- All$counts[,All$phenoData$PassScrub] # keep only cells that passed all criteria in QC
#dim(All$counts) 
All$phenoData <- All$phenoData[All$phenoData$PassScrub,]
#dim(All$phenoData) # make sure rows here matches columns in All$counts

## Now save the data again in multiple formats:
stopifnot(identical(as.character(All$featureData$FinalAnnot),rownames(All$counts)))
out <- list()
out[["counts"]] <- All$counts
out[["phenoData"]] <- All$phenoData
out[["featureData"]] <- All$featureData
saveRDS(out,file=file.path("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/", "GutEpScrubbedDFs.rds")) # this saves all of our information for counts, barcodes, and feature data as an .rds

seurat <- CreateSeuratObject(counts = All$counts, meta.data = All$phenoData) # create Seurat object of counts & pheno data
write10xCounts(x = seurat@assays$RNA@counts, path = "/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/GutEpScrubbed", version = "3") # create CellRanger-like output files of our Seurat object
saveRDS(seurat,file=file.path("/home/Jayne.Wiarda/scRNAseqEpSI/Scrublet/", "GutEpScrubbedSeurat.rds")) # this saves all of our information for counts, barcodes, and feature data as an .rds

##### Woo-hoo! We've finished all of our QC analyses! Now we move on to the fun part!

sessionInfo()

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
#[1] SeuratObject_4.0.2          Seurat_4.0.4                DropletUtils_1.12.2         SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0
#[6] Biobase_2.52.0              GenomicRanges_1.44.0        GenomeInfoDb_1.28.4         IRanges_2.26.0              S4Vectors_0.30.0           
#[11] BiocGenerics_0.38.0         MatrixGenerics_1.4.3        matrixStats_0.60.1          dplyr_1.0.7                

#loaded via a namespace (and not attached):
#[1] utf8_1.2.2                reticulate_1.21           R.utils_2.10.1            tidyselect_1.1.1          htmlwidgets_1.5.4         grid_4.1.1               
#[7] BiocParallel_1.26.2       Rtsne_0.15                ScaledMatrix_1.0.0        munsell_0.5.0             codetools_0.2-18          ica_1.0-2                
#[13] future_1.22.1             miniUI_0.1.1.1            withr_2.4.2               colorspace_2.0-2          knitr_1.34                rstudioapi_0.13          
#[19] ROCR_1.0-11               tensor_1.5                listenv_0.8.0             labeling_0.4.2            GenomeInfoDbData_1.2.6    polyclip_1.10-0          
#[25] farver_2.1.0              rhdf5_2.36.0              parallelly_1.28.1         vctrs_0.3.8               generics_0.1.0            xfun_0.26                
#[31] R6_2.5.1                  ggbeeswarm_0.6.0          graphlayouts_0.7.1        rsvd_1.0.5                locfit_1.5-9.4            miloR_1.0.0              
#[37] bitops_1.0-7              rhdf5filters_1.4.0        spatstat.utils_2.2-0      DelayedArray_0.18.0       assertthat_0.2.1          promises_1.2.0.1         
#[43] scales_1.1.1              ggraph_2.0.5              beeswarm_0.4.0            gtable_0.3.0              beachmat_2.8.1            globals_0.14.0           
#[49] goftest_1.2-2             tidygraph_1.2.0           rlang_0.4.11              splines_4.1.1             lazyeval_0.2.2            spatstat.geom_2.2-2      
#[55] broom_0.7.9               yaml_2.2.1                reshape2_1.4.4            abind_1.4-5               modelr_0.1.8              backports_1.2.1          
#[61] httpuv_1.6.3              tools_4.1.1               ggplot2_3.3.5             ellipsis_0.3.2            spatstat.core_2.3-0       RColorBrewer_1.1-2       
#[67] ggridges_0.5.3            Rcpp_1.0.7                plyr_1.8.6                sparseMatrixStats_1.4.2   zlibbioc_1.38.0           purrr_0.3.4              
#[73] RCurl_1.98-1.4            rpart_4.1-15              deldir_0.2-10             viridis_0.6.1             pbapply_1.5-0             cowplot_1.1.1            
#[79] zoo_1.8-9                 haven_2.4.3               ggrepel_0.9.1             cluster_2.1.2             fs_1.5.0                  magrittr_2.0.1           
#[85] data.table_1.14.0         scattermore_0.7           openxlsx_4.2.4            lmtest_0.9-38             reprex_2.0.1              RANN_2.6.1               
#[91] fitdistrplus_1.1-5        hms_1.1.0                 patchwork_1.1.1           mime_0.11                 evaluate_0.14             xtable_1.8-4             
#[97] rio_0.5.27                readxl_1.3.1              gridExtra_2.3             compiler_4.1.1            tibble_3.1.4              KernSmooth_2.23-20       
#[103] crayon_1.4.1              R.oo_1.24.0               htmltools_0.5.2           mgcv_1.8-36               later_1.3.0               tzdb_0.1.2               
#[109] tidyr_1.1.3               lubridate_1.7.10          DBI_1.1.1                 tweenr_1.0.2              dbplyr_2.1.1              MASS_7.3-54              
#[115] Matrix_1.3-4              car_3.0-11                readr_2.0.1               permute_0.9-5             cli_3.0.1                 R.methodsS3_1.8.1        
#[121] igraph_1.2.6              forcats_0.5.1             pkgconfig_2.0.3           foreign_0.8-81            plotly_4.9.4.1            scuttle_1.2.1            
#[127] spatstat.sparse_2.0-0     xml2_1.3.2                vipor_0.4.5               dqrng_0.3.0               XVector_0.32.0            rvest_1.0.1              
#[133] stringr_1.4.0             digest_0.6.27             sctransform_0.3.2         RcppAnnoy_0.0.19          vegan_2.5-7               spatstat.data_2.1-0      
#[139] rmarkdown_2.11            cellranger_1.1.0          leiden_0.3.9              uwot_0.1.10               edgeR_3.34.1              DelayedMatrixStats_1.14.3
#[145] curl_4.3.2                gtools_3.9.2              shiny_1.6.0               lifecycle_1.0.0           nlme_3.1-152              jsonlite_1.7.2           
#[151] Rhdf5lib_1.14.2           BiocNeighbors_1.10.0      carData_3.0-4             viridisLite_0.4.0         limma_3.48.3              fansi_0.5.0              
#[157] pillar_1.6.2              lattice_0.20-44           fastmap_1.1.0             httr_1.4.2                survival_3.2-13           glue_1.4.2               
#[163] zip_2.2.0                 png_0.1-7                 ggforce_0.3.3             stringi_1.7.4             HDF5Array_1.20.0          BiocSingular_1.8.1       
#[169] irlba_2.3.3               future.apply_1.8.1       
