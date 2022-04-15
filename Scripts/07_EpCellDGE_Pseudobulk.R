# Title: 07 - Epithelial Cell DGE Analysis & Pseudlbulk Correlations
# Date: 2021Oct06
# Author: Jayne Wiarda

library(Seurat)  
library(writexl)
library(SeuratDisk)
library(scales)
library(dplyr)
library(ggplot2)
library(readxl)

# Load in dataset:
ep <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
DefaultAssay <- 'RNA'

# Perform cluster-based DGE analysis for all clusters:
# DGE is one cluster vs average of all other epithelial cells in dataset
Idents(ep) <- ep$celltype
epDE <- FindAllMarkers(ep, # perform DGE analysis of each cluster vs average of all other cells in dataset
                       only.pos = FALSE, 
                       logfc.threshold = 0.5, # minimum logFC of 0.25
                       min.pct = 0.1, # expressed in at least 10% of all cells in the cluster or entire dataset
                       assay = "RNA") 
epDE <- subset(epDE, p_val_adj < 0.05) # make sure the adjusted p-values are still < 0.05 since some genes in DE list have p_val_adj > 0.05
pigGenes <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalAnnot <-gsub("_", "-", pigGenes$FinalAnnot) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[c('ENSID', 'FinalAnnot')]
DE <- merge(epDE, pigGenes, by.x = 'gene', by.y = 'FinalAnnot')

# Save DE gene list:
write_xlsx(x = DE, 
           path = "/home/Jayne.Wiarda/scRNAseqEpSI/DGE/EpOnly_OverallCellTypeDEGs.xlsx",
           col_names = TRUE)

# Create and store background gene list from DGE analysis (might need for later analyses):
pigGenes <- pigGenes[pigGenes$FinalAnnot %in% rownames(ep[['RNA']]@counts), ] # slim down to only genes in our dataset
write_xlsx(x = pigGenes, 
           path = "/home/Jayne.Wiarda/scRNAseqEpSI/DGE/EpOnly_BackgroundGeneList.xlsx",
           col_names = TRUE)

# Plot top DE genes:
epDE <- subset(epDE, avg_log2FC > 0) # only take genes enriched in the clusters 
#topgenes <- epDE %>% group_by(cluster) %>% top_n(30, avg_log2FC) # only plot top 50 genes per cluster, as determined by highest average log2FC values
topgenes <- epDE %>% group_by(cluster) %>% top_n(-30, p_val_adj) # find 50 genes with lowest p-values in each cluster
DoHeatmap(subset(ep, downsample = 30), # show only 30 cells from each cell type
          features = as.character(topgenes$gene), 
          assay = "RNA", 
          disp.min = -2, 
          disp.max = 2,
          draw.lines = FALSE,
          group.colors = c('burlywood4', 'goldenrod1', 'darkorange2', 'chartreuse3', 'turquoise4', 'darkmagenta')) +
  scale_fill_gradientn(colors = c('slateblue4', 'black', 'yellowgreen'))

# Perform cluster-based DGE analysis between the two enteroendocrine cell clusters:
# I prefer to use the FindAllMarkers() function over FindMarkers() since we get more output from the former, so first I'll subset to only the two clusters of interest... that way, FindAllMarkers() comparison will be equivalent to pw comparison of two clusters
ee <- subset(ep, idents = c('NEUROD1lo EE', 'NEUROD1hi EE'))

# Identify genes with non-zero expression to keep:
counts <- as.data.frame(ee[['RNA']]@counts)
keep <- rowSums(counts) > 0
keep <- rownames(counts[keep,])

# Slim down the Seurat object:
ee <- DietSeurat(ee, 
                 counts = TRUE,
                 data = TRUE,
                 scale.data = TRUE, 
                 dimreducs = NULL,
                 features = keep, # keep only genes with non-zero counts across all cells
                 assays = 'RNA') # keep only RNA assay and remove SCT and integrated
Idents(ee) <- ee$celltype
eeDE <- FindAllMarkers(ee, # perform DGE analysis of each cluster vs average of all other cells in dataset
                       only.pos = TRUE, # change to true since equivalent to pairwise comparison 
                       logfc.threshold = 0.5, # minimum logFC of 0.25
                       min.pct = 0.1, # expressed in at least 10% of all cells in the cluster or entire dataset
                       assay = "RNA") 
eeDE <- subset(eeDE, p_val_adj < 0.05) # make sure the adjusted p-values are still < 0.05 since some genes in DE list have p_val_adj > 0.05
pigGenes <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalAnnot <-gsub("_", "-", pigGenes$FinalAnnot) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[c('ENSID', 'FinalAnnot')]
DE <- merge(eeDE, pigGenes, by.x = 'gene', by.y = 'FinalAnnot')

# Save DE gene list:
write_xlsx(x = DE, 
           path = "/home/Jayne.Wiarda/scRNAseqEpSI/DGE/EpOnly_EnteroendocrinePairwiseDEGs.xlsx",
           col_names = TRUE)

# Create and store background gene list from DGE analysis (might need for later analyses):
pigGenes <- pigGenes[pigGenes$FinalAnnot %in% rownames(ee[['RNA']]@counts), ] # slim down to only genes in our dataset
write_xlsx(x = pigGenes, 
           path = "/home/Jayne.Wiarda/scRNAseqEpSI/DGE/EpOnly_EnteroendocrineSubset_BackgroundGeneList.xlsx",
           col_names = TRUE)

# Plot top DE genes:
#eeDE <- subset(eeDE, avg_log2FC > 0) # only take genes enriched in the clusters 
#topgenes <- eeDE %>% group_by(cluster) %>% top_n(50, avg_log2FC) # only plot top 50 genes per cluster, as determined by highest average log2FC values
topgenes <- eeDE %>% group_by(cluster) %>% top_n(-50, p_val_adj) # find 50 genes with lowest p-values in each cluster
DoHeatmap(ee, 
          features = as.character(topgenes$gene), 
          assay = "RNA", 
          disp.min = -2, 
          disp.max = 2,
          draw.lines = FALSE,
          group.colors = c('turquoise4', 'darkmagenta')) +
  scale_fill_gradientn(colors = c('slateblue4', 'black', 'yellowgreen'))

eeDE <- eeDE %>% group_by(cluster)
eeDE <- eeDE %>% arrange(desc(avg_log2FC), .by_group = TRUE)
DoHeatmap(ee, 
          features = eeDE$gene, # or try plotting all DE genes since not that many in this case
          assay = "RNA", 
          disp.min = -2, 
          disp.max = 2,
          draw.lines = FALSE,
          group.colors = c('turquoise4', 'darkmagenta')) +
  scale_fill_gradientn(colors = c('slateblue4', 'black', 'yellowgreen'))

# Perform pseudobulk correlation analysis on clusters:
Idents(ep) <- ep$celltype
av.exp <- AverageExpression(ep, return.seurat = TRUE) # create in-silico bulk RNA-seq dataset for each sample
avcounts <- as.matrix(av.exp@assays$RNA@data)
cor.exp <- as.data.frame(cor(avcounts, method = 'spearman'))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, levels(Idents(ep)))
cor.df$x <- factor(cor.df$x,levels = levels(Idents(ep)))
cor.df$y <- factor(cor.df$y,levels = levels(Idents(ep)))
cor.df$correlation <- round(cor.df$correlation, digits = 2)
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradientn(colours = c('yellow','orange','red', 'red3'), oob = squish, limits = c(min(cor.df$correlation), 1)) + 
  theme_classic()
ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile()+
  scale_fill_gradientn(colours = c('yellow','orange','red', 'red3'))+ 
  geom_text(aes(x, y, label = correlation), color = "black", size = 2) + theme_classic()

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
#[1] readxl_1.3.1          ggplot2_3.3.5         dplyr_1.0.7           scales_1.1.1          SeuratDisk_0.0.0.9019 writexl_1.4.0        
#[7] SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#[1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              splines_4.1.1               BiocParallel_1.26.2        
#[6] listenv_0.8.0               scattermore_0.7             GenomeInfoDb_1.28.4         digest_0.6.27               htmltools_0.5.2            
#[11] viridis_0.6.1               fansi_0.5.0                 magrittr_2.0.1              ScaledMatrix_1.0.0          tensor_1.5                 
#[16] cluster_2.1.2               ROCR_1.0-11                 limma_3.48.3                globals_0.14.0              graphlayouts_0.7.1         
#[21] matrixStats_0.60.1          spatstat.sparse_2.0-0       colorspace_2.0-2            ggrepel_0.9.1               xfun_0.26                  
#[26] crayon_1.4.1                RCurl_1.98-1.4              jsonlite_1.7.2              spatstat.data_2.1-0         ape_5.5                    
#[31] survival_3.2-13             zoo_1.8-9                   glue_1.4.2                  polyclip_1.10-0             gtable_0.3.0               
#[36] zlibbioc_1.38.0             XVector_0.32.0              leiden_0.3.9                DelayedArray_0.18.0         BiocSingular_1.8.1         
#[41] future.apply_1.8.1          SingleCellExperiment_1.14.1 BiocGenerics_0.38.0         abind_1.4-5                 edgeR_3.34.1               
#[46] DBI_1.1.1                   miniUI_0.1.1.1              Rcpp_1.0.7                  viridisLite_0.4.0           xtable_1.8-4               
#[51] reticulate_1.21             spatstat.core_2.3-0         rsvd_1.0.5                  bit_4.0.4                   htmlwidgets_1.5.4          
#[56] httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.2              ica_1.0-2                   pkgconfig_2.0.3            
#[61] farver_2.1.0                uwot_0.1.10                 dbplyr_2.1.1                deldir_0.2-10               locfit_1.5-9.4             
#[66] utf8_1.2.2                  labeling_0.4.2              tidyselect_1.1.1            rlang_0.4.11                reshape2_1.4.4             
#[71] later_1.3.0                 cellranger_1.1.0            munsell_0.5.0               tools_4.1.1                 cli_3.0.1                  
#[76] generics_0.1.0              ggridges_0.5.3              evaluate_0.14               stringr_1.4.0               fastmap_1.1.0              
#[81] yaml_2.2.1                  goftest_1.2-2               knitr_1.34                  bit64_4.0.5                 fitdistrplus_1.1-5         
#[86] tidygraph_1.2.0             purrr_0.3.4                 RANN_2.6.1                  ggraph_2.0.5                pbapply_1.5-0              
#[91] future_1.22.1               nlme_3.1-152                mime_0.11                   rstudioapi_0.13             hdf5r_1.3.4                
#[96] compiler_4.1.1              beeswarm_0.4.0              plotly_4.9.4.1              png_0.1-7                   spatstat.utils_2.2-0       
#[101] tibble_3.1.4                tweenr_1.0.2                stringi_1.7.4               lattice_0.20-45             Matrix_1.3-4               
#[106] vctrs_0.3.8                 pillar_1.6.2                lifecycle_1.0.0             spatstat.geom_2.2-2         lmtest_0.9-38              
#[111] RcppAnnoy_0.0.19            BiocNeighbors_1.10.0        data.table_1.14.0           cowplot_1.1.1               bitops_1.0-7               
#[116] irlba_2.3.3                 httpuv_1.6.3                patchwork_1.1.1             GenomicRanges_1.44.0        R6_2.5.1                   
#[121] promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3               vipor_0.4.5                 IRanges_2.26.0             
#[126] parallelly_1.28.1           codetools_0.2-18            gtools_3.9.2                MASS_7.3-54                 assertthat_0.2.1           
#[131] SummarizedExperiment_1.22.0 withr_2.4.2                 sctransform_0.3.2           S4Vectors_0.30.2            GenomeInfoDbData_1.2.6     
#[136] mgcv_1.8-37                 miloR_1.0.0                 beachmat_2.8.1              grid_4.1.1                  rpart_4.1-15               
#[141] tidyr_1.1.3                 rmarkdown_2.11              MatrixGenerics_1.4.3        Rtsne_0.15                  ggforce_0.3.3              
#[146] Biobase_2.52.0              shiny_1.6.0                 ggbeeswarm_0.6.0    