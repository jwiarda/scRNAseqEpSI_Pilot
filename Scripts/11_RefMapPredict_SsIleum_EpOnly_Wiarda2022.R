# Title: 11 - Reference Mapping and Prediction to Porcine Ileum
# Date: 2021Sept26
# Author: Jayne Wiarda

library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)
library(viridis)

# Load in an older dataset containing epithelial cells from ileum of 2 pigs
ep1 <- readRDS('/home/Jayne.Wiarda/scRNAseqIleumAtlas/Ileum/Seurat/IleumAtlasAll_FINALannot.rds')
Idents(ep1) <- ep1$cellID

# Subset down to only epithelial cells:
ep1 <- subset(ep1, idents = 'Epithelial cells')

# Convert to gene nomenclature used in our new dataset:
querygenes <- as.data.frame(rownames(ep1[['RNA']]@counts)) # extract pig gene names from dataset
colnames(querygenes) <- 'gene'
pigGenesOld <- read_excel('/home/Jayne.Wiarda/scRNAseqIleumAtlas/Ileum/QC/UnfilteredGeneInfo.xlsx')
pigGenesOld$Name <-gsub("_", "-", pigGenesOld$Name) # replace all underscores with dashes since this occurred when processing data into Seurat
pigGenesOld <- pigGenesOld[pigGenesOld$Name %in% querygenes$gene, ] # slim down to only genes in our dataset
querycounts <- ep1[['RNA']]@counts[rownames(ep1[['RNA']]@counts) %in% pigGenesOld$Name,]
pigGenesOld <- pigGenesOld %>% arrange(factor(Name, levels = rownames(querycounts))) # arrange pigGenesOld in same order as querycounts
rownames(querycounts) <- pigGenesOld$EnsemblID # change pig gene names to Ensembl IDs

pigGenes <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalAnnot <-gsub("_", "-", pigGenes$FinalAnnot) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[pigGenes$ENSID %in% rownames(querycounts), ] # slim down to only genes in our dataset
pigGenes <- pigGenes %>% arrange(factor(ENSID, levels = rownames(querycounts))) # arrange pigGenes in same order as querycounts
rownames(querycounts) <- pigGenes$FinalAnnot # change pig gene names to updated gene annotations
queryMeta <- ep1@meta.data

# Create new Seurat object
ep1 <- CreateSeuratObject(counts = querycounts, # create new Seurat object for our query dataset
                            meta.data = queryMeta)
#ep1$SampleID <- paste0('exp1_', ep1$SampleID)
rm(list= ls()[!(ls() %in% c('ep1'))]) # clear up space

# Normalization and integration:
query <- SplitObject(ep1, split.by = "SampleID") # split into the original samples that were processed for scRNA-seq
for (i in 1:length(query)) { # for each sample individually, let's normalize the data and find the 2000 most highly variable features, scale the data, and find top 50 PCs
  query[[i]] <- NormalizeData(query[[i]], 
                              verbose = TRUE, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000, 
                              assay = "RNA")
  query[[i]] <- FindVariableFeatures(query[[i]], 
                                     selection.method = "vst", 
                                     nfeatures = 2000, 
                                     verbose = TRUE)
  query[[i]] <- ScaleData(query[[i]]) # scale the data
  query[[i]] <- RunPCA(query[[i]], # calculate 30 PCs
                       npcs = 40, 
                       verbose = TRUE)
}

# Now save the query data:
saveRDS(query, '/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/Wiarda_IleumAtlas_query.rds')

# Load in our data from newer dataset where epithelial cells are further annotated:
ref <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
DefaultAssay(ref) <- 'integrated'
Idents(ref) <- ref$celltype
#ep2$SampleID <- paste0('exp2_', ep2$orig.ident)

## Perform label transfer and mapping:
CellTypePredictions <- list()
for(i in 1:length(query)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query[[i]],
    #reduction = "cca", # opted to use cca since the method is recommended for cross-species mapping 
    reduction = 'pcaproject', # not using cca unless cross-species comparison is being performed
    dims = 1:30, 
    normalization.method = "SCT",
    k.filter = 46) # reduce k.filter from 200 to 46 since our smallest sample only has 46 cells
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$celltype), 
                              dims = 1:30,
                              weight.reduction = "pcaproject",
                              k.weight = 25)
  CellTypePredictions[[i]] <- predictions
} 

CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)

# Save the mapping & prediction results:
CellTypePredictions$CellBarcodes <- rownames(CellTypePredictions)
write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/WiardaIleumAtals_EpOnly_CellTypePredictions.xlsx')

# Incorporate cell prediction & mapping scores into original Seurat object of query data:
ep <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/PorcineIleumAtlasGut1Experiments_EpOnly.h5Seurat')
ep <- AddMetaData(object = ep, 
                  metadata = c(CellTypePredictions))

predict <- colnames(ep@meta.data %>% select(starts_with("prediction.score."))) # extract 
names <- sub(".*prediction.score. ", "", predict)   
FeaturePlot(ep,
            features = c(predict),
            reduction = 'tsne',
            ncol = 5) & 
  scale_color_gradientn( colours = c('grey90', 'darkslateblue'),  limits = c(0, 1), oob = squish) & 
  NoAxes() & NoLegend() 
FeaturePlot(ep,
            features = c(predict),
            reduction = 'tsne',
            ncol = 5) & 
  scale_color_viridis() & 
  NoAxes() & NoLegend() 
names

DimPlot(ep, 
        reduction = 'tsne',
        group.by = 'predicted.id', cols = c('darkorange2',
                                            'burlywood4',
                                            'goldenrod1',
                                            'chartreuse3',
                                            'darkmagenta',
                                            'turquoise4'))

sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.4 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#[6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] scales_1.1.1          writexl_1.4.0         readxl_1.3.1          dplyr_1.0.7           SeuratDisk_0.0.0.9019 Seurat_4.0.4         
#[7] ggplot2_3.3.5         viridis_0.6.2         viridisLite_0.4.0     SeuratObject_4.0.2   

#loaded via a namespace (and not attached):
#  [1] nlme_3.1-157          matrixStats_0.60.1    spatstat.sparse_2.0-0 bit64_4.0.5           RcppAnnoy_0.0.19      RColorBrewer_1.1-2   
#[7] httr_1.4.2            sctransform_0.3.2     tools_4.1.3           utf8_1.2.2            R6_2.5.1              irlba_2.3.3          
#[13] rpart_4.1.16          KernSmooth_2.23-20    uwot_0.1.10           mgcv_1.8-40           DBI_1.1.1             lazyeval_0.2.2       
#[19] colorspace_2.0-2      withr_2.4.2           tidyselect_1.1.1      gridExtra_2.3         bit_4.0.4             compiler_4.1.3       
#[25] cli_3.0.1             hdf5r_1.3.4           plotly_4.9.4.1        labeling_0.4.2        lmtest_0.9-38         spatstat.data_2.1-0  
#[31] ggridges_0.5.3        pbapply_1.5-0         goftest_1.2-2         stringr_1.4.0         digest_0.6.27         spatstat.utils_2.2-0 
#[37] pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.28.1     fastmap_1.1.0         htmlwidgets_1.5.4     rlang_0.4.11         
#[43] shiny_1.6.0           farver_2.1.0          generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2            
#[49] magrittr_2.0.1        patchwork_1.1.1       Matrix_1.4-1          Rcpp_1.0.7            munsell_0.5.0         fansi_0.5.0          
#[55] abind_1.4-5           reticulate_1.21       lifecycle_1.0.0       stringi_1.7.4         MASS_7.3-56           Rtsne_0.15           
#[61] plyr_1.8.6            grid_4.1.3            parallel_4.1.3        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1        
#[67] crayon_1.4.1          deldir_0.2-10         miniUI_0.1.1.1        lattice_0.20-45       cowplot_1.1.1         splines_4.1.3        
#[73] tensor_1.5            pillar_1.6.2          igraph_1.2.6          spatstat.geom_2.2-2   future.apply_1.8.1    reshape2_1.4.4       
#[79] codetools_0.2-18      leiden_0.3.9          glue_1.4.2            data.table_1.14.0     png_0.1-7             vctrs_0.3.8          
#[85] httpuv_1.6.3          cellranger_1.1.0      gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0  
#[91] polyclip_1.10-0       tidyr_1.1.3           scattermore_0.7       future_1.22.1         assertthat_0.2.1      mime_0.11            
#[97] xtable_1.8-4          later_1.3.0           survival_3.3-1        tibble_3.1.4          cluster_2.1.3         globals_0.14.0       
#[103] fitdistrplus_1.1-5    ellipsis_0.3.2        ROCR_1.0-11          