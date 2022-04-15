# Title: 15 - Plots of Enriched Genes and Processes
# Date: 2022April15
# Author: Jayne Wiarda

library(Seurat)
library(SeuratDisk)
library(readxl)
library(ggplot2)

# Load in dataset:
ep <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
DefaultAssay <- 'RNA'
Idents(ep) <- ep$celltype

# Create dot plot of selected DE genes
## Genes selected from top 25 genes based on highest logFC values
DotPlot(ep, 
        features = c('RPL5', 'RPS6', 'EEF1B2', 'OLFM4', 'PIGR', 'LYZ', 
                     'FABP2', 'FABP1', 'CLCA4', 'SLC5A1', 'SI', 'ACE2',
                     'GUCA2A', 'GUCA2B', 'BEST4', 'CFTR', 'NOTCH2', 'OTOP2', 
                     'TFF3', 'REG4', 'CLCA1', 'SPINK4', 'MUC2', 'CXCL8', 
                      'PYY', 'GAST', 'SST', 'CCK','TTR', 'NTS', 
                     'NEUROD1', 'CHGA', 'CHGB', 'KRT7', 'SCT', 'PENK'),
        cols = c('gold', 'darkgreen')) + RotatedAxis()

## Only hormone-encoding genes in enteroendocrine cells
ee <- subset(ep, idents = c('NEUROD1lo EE', 'NEUROD1hi EE'))
ee$combo <- paste(ee$celltype, ee$region, sep = ' - ')
Idents(ee) <- ee$combo
levels(ee)
FeaturePlot(ee, 
            features = c('SCT', 'PENK', 'SST', 'CCK', 'GAST', 'TTR',
                         'GHRL', 'MLN', 'GIP', 'PYY', 'NTS', 'GCG'),
            cols = c('grey80', 'darkgreen'),
            pt.size = .5) & NoAxes() & NoLegend() 
levels(ee) <- c('NEUROD1lo EE - Duodenum', 'NEUROD1lo EE - Jejunum', 'NEUROD1lo EE - Ileum',
                'NEUROD1hi EE - Duodenum', 'NEUROD1hi EE - Jejunum', 'NEUROD1hi EE - Ileum')
DotPlot(ee, 
        features = c('SCT', 'PENK', 'SST', 'CCK', 'GAST', 'TTR',
                     'GHRL', 'MLN', 'GIP', 'PYY', 'NTS', 'GCG'),
        cols = c('yellow', 'darkgreen')) + RotatedAxis()

# Also re-make t-SNE of only EE cells
DimPlot(ee,
        group.by = 'celltype',
        reduction = 'tsne',
        cols = c('turquoise4', 'darkmagenta'))
DimPlot(ee,
        group.by = 'region',
        reduction = 'tsne')

# Create dot plot of selected GO processes for all epithelial cell types
GO <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/GO/Ep2_GO_AllCells.xlsx')
head(GO) # need to fourth row column names and then remove first four rows
colnames(GO) <- GO[4,]
GO <- GO[-c(1:4),]
colnames(GO)
GO$FDR <- as.numeric(GO$FDR)
GO$fold_enrichment <- as.numeric(GO$fold_enrichment)
names(GO) <- gsub(" ", "_", names(GO)) # remove spaces in names
terms <- c('positive regulation of cell migration (GO:0030335)',
           'positive regulation of cell motility (GO:2000147)', 
           'positive regulation of response to stimulus (GO:0048584)',
           'cell differentiation (GO:0030154)',
           'cytoplasmic translation (GO:0002181)',
           'ribosomal small subunit assembly (GO:0000028)',
           'ATP synthesis coupled proton transport (GO:0015986)',
           'cell division (GO:0051301)',
           'biosynthetic process (GO:0009058)',
           'positive regulation of cell cycle (GO:0045787)',
           'mitotic cell cycle (GO:0000278)',
           'bile acid secretion (GO:0032782)',
           'modified amino acid transport (GO:0072337)',
           'fatty acid beta-oxidation (GO:0006635)',
           'lipid oxidation (GO:0034440)',
           'import across plasma membrane (GO:0098739)',
           'sodium ion transport (GO:0006814)',
           'posttranslational protein targeting to membrane, translocation (GO:0031204)',
           'nucleotide-sugar biosynthetic process (GO:0009226)',
           'protein localization to endoplasmic reticulum (GO:0070972)',
           'Golgi vesicle transport (GO:0048193)',
           'glycosylation (GO:0070085)',
           'glycoprotein biosynthetic process (GO:0009101)',
           'regulation of secretion by cell (GO:1903530)',
           'regulation of hormone levels (GO:0010817)')
GO <- GO[GO$GO_biological_process_complete %in% terms, ]
GO$cell_type <- factor(GO$cell_type, 
                       levels=c('crypt', 'enterocyte', 'BEST4 enterocyte',
                                'goblet', 'NEUROD1lo EE', 'NEUROD1hi EE'))
GO$GO_biological_process_complete <- factor(GO$GO_biological_process_complete,
                                            levels = c('cytoplasmic translation (GO:0002181)',
                                                       'ribosomal small subunit assembly (GO:0000028)',
                                                       'ATP synthesis coupled proton transport (GO:0015986)',
                                                       'cell division (GO:0051301)',
                                                       'biosynthetic process (GO:0009058)',
                                                       'positive regulation of cell cycle (GO:0045787)',
                                                       'mitotic cell cycle (GO:0000278)',
                                                       'bile acid secretion (GO:0032782)',
                                                       'modified amino acid transport (GO:0072337)',
                                                       'fatty acid beta-oxidation (GO:0006635)',
                                                       'lipid oxidation (GO:0034440)',
                                                       'import across plasma membrane (GO:0098739)',
                                                       'sodium ion transport (GO:0006814)',
                                                       'positive regulation of cell migration (GO:0030335)',
                                                       'positive regulation of cell motility (GO:2000147)', 
                                                       'positive regulation of response to stimulus (GO:0048584)',
                                                       'cell differentiation (GO:0030154)',
                                                       'posttranslational protein targeting to membrane, translocation (GO:0031204)',
                                                       'nucleotide-sugar biosynthetic process (GO:0009226)',
                                                       'protein localization to endoplasmic reticulum (GO:0070972)',
                                                       'Golgi vesicle transport (GO:0048193)',
                                                       'glycosylation (GO:0070085)',
                                                       'glycoprotein biosynthetic process (GO:0009101)',
                                                       'regulation of secretion by cell (GO:1903530)',
                                                       'regulation of hormone levels (GO:0010817)'))
ggplot(GO) +
  geom_point(aes(x = GO_biological_process_complete, 
                 y = cell_type,
                 size = fold_enrichment,
                 color = FDR)) +
  theme_bw() +
  scale_color_gradient(low = 'darkgreen', high = 'gold') +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

# Create dot plot of selected GO processes for only EE cells
GO <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/GO/Ep2_GO_EE.xlsx')
head(GO) # need to fourth row column names and then remove first four rows
colnames(GO) <- GO[4,]
GO <- GO[-c(1:4),]
colnames(GO)
GO$FDR <- as.numeric(GO$FDR)
GO$fold_enrichment <- as.numeric(GO$fold_enrichment)
names(GO) <- gsub(" ", "_", names(GO)) # remove spaces in names
terms <- c('mitochondrial ATP synthesis coupled electron transport (GO:0042775)',
           'aerobic respiration (GO:0009060)',
           'ribosome assembly (GO:0042255)',
           'cytoplasmic translation (GO:0002181)',
           'peptide biosynthetic process (GO:0043043)')
GO <- GO[GO$GO_biological_process_complete %in% terms, ]
GO$cell_type <- factor(GO$cell_type, 
                       levels=c('NEUROD1lo EE', 'NEUROD1hi EE'))
GO$GO_biological_process_complete <- factor(GO$GO_biological_process_complete,
                                            levels = c('mitochondrial ATP synthesis coupled electron transport (GO:0042775)',
                                                       'aerobic respiration (GO:0009060)',
                                                       'ribosome assembly (GO:0042255)',
                                                       'cytoplasmic translation (GO:0002181)',
                                                       'peptide biosynthetic process (GO:0043043)'))
ggplot(GO) +
  geom_point(aes(x = GO_biological_process_complete, 
                 y = cell_type,
                 size = fold_enrichment,
                 color = FDR)) +
  theme_bw() +
  scale_color_gradient(low = 'darkgreen', high = 'gold') +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))

sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggplot2_3.3.5         readxl_1.3.1          SeuratDisk_0.0.0.9019 SeuratObject_4.0.2    Seurat_4.0.4         

#loaded via a namespace (and not attached):
#  [1] utf8_1.2.2                  reticulate_1.21             R.utils_2.10.1              tidyselect_1.1.1           
#[5] RSQLite_2.2.8               AnnotationDbi_1.54.1        htmlwidgets_1.5.4           grid_4.1.3                 
#[9] Rtsne_0.15                  munsell_0.5.0               codetools_0.2-18            ica_1.0-2                  
#[13] future_1.22.1               miniUI_0.1.1.1              withr_2.4.2                 colorspace_2.0-2           
#[17] Biobase_2.52.0              rstudioapi_0.13             stats4_4.1.3                ROCR_1.0-11                
#[21] tensor_1.5                  listenv_0.8.0               labeling_0.4.2              MatrixGenerics_1.4.3       
#[25] GenomeInfoDbData_1.2.6      polyclip_1.10-0             bit64_4.0.5                 farver_2.1.0               
#[29] parallelly_1.28.1           vctrs_0.3.8                 generics_0.1.0              R6_2.5.1                   
#[33] GenomeInfoDb_1.28.4         hdf5r_1.3.4                 bitops_1.0-7                spatstat.utils_2.2-0       
#[37] cachem_1.0.6                DelayedArray_0.18.0         assertthat_0.2.1            promises_1.2.0.1           
#[41] scales_1.1.1                gtable_0.3.0                globals_0.14.0              goftest_1.2-2              
#[45] rlang_0.4.11                splines_4.1.3               lazyeval_0.2.2              spatstat.geom_2.2-2        
#[49] broom_0.7.9                 reshape2_1.4.4              abind_1.4-5                 modelr_0.1.8               
#[53] backports_1.2.1             httpuv_1.6.3                tools_4.1.3                 ellipsis_0.3.2             
#[57] spatstat.core_2.3-0         RColorBrewer_1.1-2          BiocGenerics_0.38.0         ggridges_0.5.3             
#[61] Rcpp_1.0.7                  plyr_1.8.6                  zlibbioc_1.38.0             purrr_0.3.4                
#[65] RCurl_1.98-1.4              rpart_4.1.16                deldir_0.2-10               pbapply_1.5-0              
#[69] cowplot_1.1.1               S4Vectors_0.30.2            zoo_1.8-9                   SummarizedExperiment_1.22.0
#[73] haven_2.4.3                 ggrepel_0.9.1               cluster_2.1.3               fs_1.5.0                   
#[77] magrittr_2.0.1              data.table_1.14.0           scattermore_0.7             openxlsx_4.2.4             
#[81] lmtest_0.9-38               reprex_2.0.1                RANN_2.6.1                  fitdistrplus_1.1-5         
#[85] matrixStats_0.60.1          hms_1.1.0                   patchwork_1.1.1             mime_0.11                  
#[89] xtable_1.8-4                XML_3.99-0.7                rio_0.5.27                  IRanges_2.26.0             
#[93] gridExtra_2.3               compiler_4.1.3              tibble_3.1.4                KernSmooth_2.23-20         
#[97] crayon_1.4.1                R.oo_1.24.0                 htmltools_0.5.2             mgcv_1.8-40                
#[101] later_1.3.0                 tzdb_0.1.2                  tidyr_1.1.3                 lubridate_1.7.10           
#[105] DBI_1.1.1                   dbplyr_2.1.1                MASS_7.3-56                 Matrix_1.4-1               
#[109] car_3.0-11                  readr_2.0.1                 cli_3.0.1                   R.methodsS3_1.8.1          
#[113] parallel_4.1.3              igraph_1.2.6                GenomicRanges_1.44.0        forcats_0.5.1              
#[117] pkgconfig_2.0.3             foreign_0.8-82              plotly_4.9.4.1              spatstat.sparse_2.0-0      
#[121] xml2_1.3.2                  annotate_1.70.0             XVector_0.32.0              rvest_1.0.1                
#[125] stringr_1.4.0               digest_0.6.27               sctransform_0.3.2           RcppAnnoy_0.0.19           
#[129] graph_1.70.0                spatstat.data_2.1-0         Biostrings_2.60.2           cellranger_1.1.0           
#[133] leiden_0.3.9                uwot_0.1.10                 GSEABase_1.54.0             curl_4.3.2                 
#[137] shiny_1.6.0                 lifecycle_1.0.0             nlme_3.1-157                jsonlite_1.7.2             
#[141] carData_3.0-4               viridisLite_0.4.0           fansi_0.5.0                 pillar_1.6.2               
#[145] lattice_0.20-45             KEGGREST_1.32.0             fastmap_1.1.0               httr_1.4.2                 
#[149] survival_3.3-1              glue_1.4.2                  zip_2.2.0                   png_0.1-7                  
#[153] bit_4.0.4                   stringi_1.7.4               blob_1.2.2                  memoise_2.0.0              
#[157] dplyr_1.0.7                 irlba_2.3.3                 future.apply_1.8.1   