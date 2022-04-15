# Title: 12 - Gene Set Enrichment Analysis w/ Porcine Ileum Dataset
# Date: 2022April05
# Author: Jayne Wiarda

library(Seurat)
library(SeuratDisk)
library(AUCell)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(writexl)
library(readxl)
library(viridis)

# Create gene signatures from DE genes:
## Create gene sets:
DE <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/DGE/EpOnly_OverallCellTypeDEGs.xlsx')
DE <- subset(DE, avg_log2FC > 0) # keep only genes with positive fold changes
crypt <- subset(DE, cluster == 'Crypt')
crypt <- crypt$gene
ent <- subset(DE, cluster == 'Enterocyte')
ent <- ent$gene
best4 <- subset(DE, cluster == 'BEST4 enterocyte')
best4 <- best4$gene
gob <- subset(DE, cluster == 'Goblet')
gob <- gob$gene
EElo <- subset(DE, cluster == 'NEUROD1lo EE')
EElo <- EElo$gene
EEhi <- subset(DE, cluster == 'NEUROD1hi EE')
EEhi <- EEhi$gene

## Create list of all gene sets:
geneSets <- list(Crypt = crypt, Enterocyte = ent,
                 BEST4_enterocyte = best4, Goblet = gob,
                 NEUROD1lo_EE = EElo, NEUROD1hi_EE = EEhi)

# Perform gene set enrichment analysis:

##Create experimental matrix from raw counts data:

ep.integrated <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/PorcineIleumAtlasGut1Experiments_EpOnly.h5Seurat')
exprMatrix <- as.matrix(ep.integrated[["RNA"]]@counts) 

#Calculate gene rankings within each cell…this can be used to calculate enrichment scores for any gene sets:
  
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE) 

#Calculate AUC scores for each cell with each gene set… we use the default top 5% of expressed genes to calculate our AUC score:
  
#cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # calculate cell signature AUC score for each gene set in each cell
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings))) # to change top % of expressed genes used to caluclate AUC to 5%

#Collate results:
  
AUCs <- as.data.frame(getAUC(cells_AUC))
AUCs <- t(AUCs)
AUCs <- as.data.frame(AUCs)
AUCs <- mutate_all(AUCs, function(x) as.numeric(as.character(x)))
head(AUCs)

# Save the AUC scores:
AUCs$CellBarcodes <- rownames(AUCs)
write_xlsx(AUCs, '/home/Jayne.Wiarda/scRNAseqEpSI/AUC/WiardaIleumAtals_EpOnly_AUCscores.xlsx')

# Overlay onto tSNE:

# Crypt cells:
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$Crypt)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# Enterocytes
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$Enterocyte)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# BEST4 enterocytes
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$BEST4_enterocyte)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# Goblet cells
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$Goblet)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# NEUROD1lo EE
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$NEUROD1lo_EE)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# NEUROD1hi EE
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$NEUROD1hi_EE)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# Compare back to cell type prediction gradients:
CellTypePredictions <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/WiardaIleumAtals_EpOnly_CellTypePredictions.xlsx')
colnames(CellTypePredictions)

# Crypt cells:
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, CellTypePredictions$prediction.score.Crypt)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# Enterocytes
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, CellTypePredictions$prediction.score.Enterocyte)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# BEST4 enterocytes
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, CellTypePredictions$prediction.score.BEST4.enterocyte)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# Goblet cells
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, CellTypePredictions$prediction.score.Goblet)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# NEUROD1lo EE
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, CellTypePredictions$prediction.score.NEUROD1lo.EE)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

# NEUROD1hi EE
coords <- Embeddings(ep.integrated[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, CellTypePredictions$prediction.score.NEUROD1hi.EE)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis') +
  theme_classic()

sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] viridis_0.6.2         viridisLite_0.4.0     SeuratDisk_0.0.0.9019 readxl_1.3.1          writexl_1.4.0        
#[6] scales_1.1.1          tidyr_1.1.3           ggplot2_3.3.5         dplyr_1.0.7           AUCell_1.14.0        
#[11] SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#  [1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              GSEABase_1.54.0            
#[5] splines_4.1.3               listenv_0.8.0               scattermore_0.7             GenomeInfoDb_1.28.4        
#[9] digest_0.6.27               htmltools_0.5.2             fansi_0.5.0                 magrittr_2.0.1             
#[13] memoise_2.0.0               tensor_1.5                  cluster_2.1.3               ROCR_1.0-11                
#[17] globals_0.14.0              Biostrings_2.60.2           annotate_1.70.0             matrixStats_0.60.1         
#[21] R.utils_2.10.1              spatstat.sparse_2.0-0       colorspace_2.0-2            blob_1.2.2                 
#[25] ggrepel_0.9.1               crayon_1.4.1                RCurl_1.98-1.4              jsonlite_1.7.2             
#[29] graph_1.70.0                spatstat.data_2.1-0         survival_3.3-1              zoo_1.8-9                  
#[33] glue_1.4.2                  polyclip_1.10-0             gtable_0.3.0                zlibbioc_1.38.0            
#[37] XVector_0.32.0              leiden_0.3.9                DelayedArray_0.18.0         future.apply_1.8.1         
#[41] BiocGenerics_0.38.0         abind_1.4-5                 DBI_1.1.1                   miniUI_0.1.1.1             
#[45] Rcpp_1.0.7                  xtable_1.8-4                reticulate_1.21             spatstat.core_2.3-0        
#[49] bit_4.0.4                   stats4_4.1.3                htmlwidgets_1.5.4           httr_1.4.2                 
#[53] RColorBrewer_1.1-2          ellipsis_0.3.2              ica_1.0-2                   pkgconfig_2.0.3            
#[57] XML_3.99-0.7                R.methodsS3_1.8.1           farver_2.1.0                uwot_0.1.10                
#[61] deldir_0.2-10               utf8_1.2.2                  tidyselect_1.1.1            labeling_0.4.2             
#[65] rlang_0.4.11                reshape2_1.4.4              later_1.3.0                 AnnotationDbi_1.54.1       
#[69] munsell_0.5.0               cellranger_1.1.0            tools_4.1.3                 cachem_1.0.6               
#[73] cli_3.0.1                   generics_0.1.0              RSQLite_2.2.8               ggridges_0.5.3             
#[77] stringr_1.4.0               fastmap_1.1.0               goftest_1.2-2               bit64_4.0.5                
#[81] fitdistrplus_1.1-5          purrr_0.3.4                 RANN_2.6.1                  KEGGREST_1.32.0            
#[85] pbapply_1.5-0               future_1.22.1               nlme_3.1-157                mime_0.11                  
#[89] R.oo_1.24.0                 hdf5r_1.3.4                 compiler_4.1.3              rstudioapi_0.13            
#[93] plotly_4.9.4.1              png_0.1-7                   spatstat.utils_2.2-0        tibble_3.1.4               
#[97] stringi_1.7.4               lattice_0.20-45             Matrix_1.4-1                vctrs_0.3.8                
#[101] pillar_1.6.2                lifecycle_1.0.0             spatstat.geom_2.2-2         lmtest_0.9-38              
#[105] RcppAnnoy_0.0.19            data.table_1.14.0           cowplot_1.1.1               bitops_1.0-7               
#[109] irlba_2.3.3                 GenomicRanges_1.44.0        httpuv_1.6.3                patchwork_1.1.1            
#[113] R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3              
#[117] IRanges_2.26.0              parallelly_1.28.1           codetools_0.2-18            MASS_7.3-56                
#[121] assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.4.2                 sctransform_0.3.2          
#[125] S4Vectors_0.30.2            GenomeInfoDbData_1.2.6      mgcv_1.8-40                 parallel_4.1.3             
#[129] grid_4.1.3                  rpart_4.1.16                MatrixGenerics_1.4.3        Rtsne_0.15                 
#[133] Biobase_2.52.0              shiny_1.6.0    