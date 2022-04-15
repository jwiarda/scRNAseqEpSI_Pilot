# Title: 13 - Gene Set Enrichment Analysis w/ Murine EE Cell Signatures
# Date: 2022April07
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

# Create gene signatures from gene lists:
## Read in gene list
DE <- read.csv('/home/Jayne.Wiarda/scRNAseqEpSI/RefData/Gehart2019_GeneListsEECs/1-s2.0-S009286741831643X-mmc3.csv')

## Load in ortho genes
orthoGenes <- read.delim("/home/Jayne.Wiarda/scRNAseqEpSI/PigToMouse_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Mouse.homology.type == 'ortholog_one2one') # subset to only one to one orthologs
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Mouse.gene.stable.ID, .keep_all = TRUE)  # retain only unique combinations of pig & human Ensembl IDs, ignoring transcript IDs

## Load in pig gene symbol key
pigGenes <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalAnnot <-gsub("_", "-", pigGenes$FinalAnnot) # replace all underscores with dashes since this occurred when processing data in a previous step

## Create gene sets
L <- DE$L.cells
L <- intersect(L, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
L <- orthoGenes[orthoGenes$Mouse.gene.name %in% L, ]
L <- pigGenes[pigGenes$ENSID %in% L$Gene.stable.ID, ] # slim down to only genes in our dataset

I <- DE$I.cells
I <- intersect(I, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
I <- orthoGenes[orthoGenes$Mouse.gene.name %in% I, ]
I <- pigGenes[pigGenes$ENSID %in% I$Gene.stable.ID, ] # slim down to only genes in our dataset

X <- DE$X.cells
X <- intersect(X, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
X <- orthoGenes[orthoGenes$Mouse.gene.name %in% X, ]
X <- pigGenes[pigGenes$ENSID %in% X$Gene.stable.ID, ] # slim down to only genes in our dataset

K <- DE$K.cells
K <- intersect(K, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
K <- orthoGenes[orthoGenes$Mouse.gene.name %in% K, ]
K <- pigGenes[pigGenes$ENSID %in% K$Gene.stable.ID, ] # slim down to only genes in our dataset

Delta <- DE$Delta.cells
Delta <- intersect(Delta, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
Delta <- orthoGenes[orthoGenes$Mouse.gene.name %in% Delta, ]
Delta <- pigGenes[pigGenes$ENSID %in% Delta$Gene.stable.ID, ] # slim down to only genes in our dataset

N <- DE$N.cells
N <- intersect(N, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
N <- orthoGenes[orthoGenes$Mouse.gene.name %in% N, ]
N <- pigGenes[pigGenes$ENSID %in% N$Gene.stable.ID, ] # slim down to only genes in our dataset

ECearly <- DE$EC.cells..early.
ECearly <- intersect(ECearly, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
ECearly <- orthoGenes[orthoGenes$Mouse.gene.name %in% ECearly, ]
ECearly <- pigGenes[pigGenes$ENSID %in% ECearly$Gene.stable.ID, ] # slim down to only genes in our dataset

EClate <- DE$EC.cells..late.
EClate <- intersect(EClate, orthoGenes$Mouse.gene.name) # find which gene names from reference are also one-to-one orthologs
EClate <- orthoGenes[orthoGenes$Mouse.gene.name %in% EClate, ]
EClate <- pigGenes[pigGenes$ENSID %in% EClate$Gene.stable.ID, ] # slim down to only genes in our dataset

## Create list of all gene sets:
geneSets <- list(L_cell = L$FinalAnnot, I_cell = I$FinalAnnot,
                 X_cell = X$FinalAnnot, K_cell = K$FinalAnnot,
                 Delta_cell = Delta$FinalAnnot, N_cell = N$FinalAnnot,
                 EC_early = ECearly$FinalAnnot, EC_late = EClate$FinalAnnot)

# Perform gene set enrichment analysis:

##Create experimental matrix from raw counts data:

ep.integrated <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')

### Use only EE cells and genes expressed by EE cells:
Idents(ep.integrated) <- ep.integrated$celltype
ee <- subset(ep.integrated, idents = c('NEUROD1lo EE', 'NEUROD1hi EE')) # subset to only EE cells 
rm(ep.integrated)

### Create matrix:
exprMatrix <- as.matrix(ee[["RNA"]]@counts) 
exprMatrix <- exprMatrix[rowSums(exprMatrix[])>0,] # remove all genes with sum zero counts

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
write_xlsx(AUCs, '/home/Jayne.Wiarda/scRNAseqEpSI/AUC/EpSI_EEonly_AUCscores.xlsx')

# Overlay onto tSNE:

# L cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$L_cell)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic() 

# I cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$I_cell)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# X cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$X_cell)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# K cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$K_cell)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# Delta cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$Delta_cell)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# N cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$N_cell)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# EC early cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$EC_early)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# EC late cells:
coords <- Embeddings(ee[["tsne"]])
coords <- as.data.frame(coords)
coords <- cbind(coords, AUCs$EC_late)
ggplot(coords, aes(x = tSNE_1, y = tSNE_2, color = coords[,3])) +
  geom_point() +
  scale_color_viridis(option = 'viridis', limits = c(0, .25), oob = squish) +
  theme_classic()

# Compare AUC scores for different modules across different annotations:
ee <- AddMetaData(object = ee, 
                  metadata = c(AUCs))

VlnPlot(ee, 
        features = colnames(AUCs[1:8]),
        cols = c('turquoise4', 'darkmagenta'),
        pt.size = 0.01,
        group.by = 'celltype',
        ncol = 4) & 
  ylim(c(0,.4)) &
  geom_boxplot(width = 0.15, fill = 'white', col = 'black', lwd = 0.5, outlier.shape = NA, coef = 0) &
  stat_summary(fun.y = mean, geom='point', size = 2, colour = "red") # red dot shows mean; box shows IQR

# Run statistics on AUC scores
meta <- ee@meta.data
meta <- meta[, c('celltype', c(colnames(AUCs[1:8])))]

## Mann-Whitney tehttp://arsiaame0fsep11.usda.net:8787/graphics/31222504-c915-4f75-8dd2-6594b6d11b63.pngsts
wilcox.test(L_cell ~ celltype, data=meta) 
wilcox.test(I_cell ~ celltype, data=meta) 
wilcox.test(X_cell ~ celltype, data=meta) 
wilcox.test(K_cell ~ celltype, data=meta) 
wilcox.test(Delta_cell ~ celltype, data=meta) 
wilcox.test(N_cell ~ celltype, data=meta) 
wilcox.test(EC_early ~ celltype, data=meta) 
wilcox.test(EC_late ~ celltype, data=meta) 

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
