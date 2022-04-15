# Title: 09 - Reference Mapping and Prediction to Human Small Intestine
# Date: 2021Sept26
# Author: Jayne Wiarda

library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)

# Load in ortho genes
orthoGenes <- read.delim("/home/Jayne.Wiarda/scRNAseqEpSI/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs

# Load in reference data
Convert( # convert downloaded .h5ad formatted dataset into an .h5seurat object. This will be deposited into the same directory as the .h5ad file
  "/home/Jayne.Wiarda/scRNAseqEpSI/RefData/Elmentaite2021_epi_raw_counts02.h5ad",
  dest = "h5seurat",
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE)
ref <- LoadH5Seurat("/home/Jayne.Wiarda/scRNAseqEpSI/RefData/Elmentaite2021_epi_raw_counts02.h5seurat") # load in .h5seurat file 

# Slim down reference data
#ref # let's see what we have in our Seurat object...
Idents(ref) <- ref$Region
#levels(Idents(ref))
ref <- subset(ref, idents = "SmallInt") # subset to take only cells from small intestine (omit large intestine, REC, APD, lymph node)
Idents(ref) <- ref$Diagnosis
#levels(Idents(ref))
ref <- subset(ref, idents = c("Pediatric healthy", "Healthy adult")) # subset to take only cells from healthy adults & pediatrics (omit Pediatric Crohn's and fetal samples)
Idents(ref) <- ref$Fraction
#levels(Idents(ref))
ref <- subset(ref, idents = c('SC-45N', 'SC-EPCAMP', 'SC')) # choosing to omit SC-45P & SC-EPCAMN fractions from our data since these should not contain large amounts of epithelial cells and many of these samples have too few cells to allow robust integration
#ref # see we have reduced cell numbers
#table(ref$sample.name) # see how many cells in each sample...our smallest sample is 38 cells, which we will need to know for setting some integration parameters below

# Filter to only one-to-one gene orthologs:
refgenes <- rownames(ref[['RNA']]@counts) # extract gene names from reference dataset
reforthos <- intersect(refgenes, orthoGenes$Human.gene.name) # find which gene names from reference are also one-to-one orthologs
#length(reforthos) # how many genes are orthologs?

# Create new Seurat object:
refcounts <- ref[['RNA']]@counts[rownames(ref[['RNA']]@counts) %in% reforthos,] # make count matrix from referemce, only taking counts from one-to-one ortholog genes
refMeta <- ref@meta.data # extract all the meta data from reference
ref <- CreateSeuratObject( # now create new Seurat object with only the filtered cells, orthologous genes, and meta data for filtered cells
  counts = refcounts, 
  meta.data = refMeta)

# Normalize & integrate data:
ref.list <- SplitObject(ref, split.by = "sample.name") # split into the original samples that were processed for scRNA-seq
for (i in 1:length(ref.list)) { # for each sample individually, let's normalize the data and find the 2000 most highly variable features
  ref.list[[i]] <- NormalizeData(ref.list[[i]], 
                                 verbose = TRUE, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000, 
                                 assay = "RNA")
  ref.list[[i]] <- FindVariableFeatures(ref.list[[i]], 
                                        selection.method = "vst", 
                                        nfeatures = 2000, 
                                        verbose = TRUE)
}
ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      dims = 1:30) # find integration anchors between samples based on variable features for each sample
ref <- IntegrateData(anchorset = ref.anchors, 
                     k.weight = 38, # reduce k weight from 100 to 38 since our smallest sample has 38 cells
                     dims = 1:30) # integrate the data together based on integration anchors found with default parameters
ref <- ScaleData(ref, 
                 verbose = TRUE, 
                 assay = 'integrated') # scale the genes in the integrated assay

# Calculate principle components:
ref <- RunPCA(ref, # calculate first 100 PCs
              npcs = 100, 
              verbose = TRUE)
ElbowPlot(ref,
          ndims = 100) # look at this plot to find the 'elbow' for significant PCs... use this number of PCs for creating UMAP, tSNE, & cell neighbors & clustering
pct <- ref[["pca"]]@stdev / sum(ref[["pca"]]@stdev) * 100 # find standard deviation for each PC
cumu <- cumsum(pct) # find cumulative percentages for PCs
co1 <- which(cumu > 90 & pct < 5)[1] # find PC representing cumulative percent >90% and less than 5% associated with the single PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # find last PC where change in percent variation is more than 0.1%
pcs <- min(co1, co2) # find the minimum PC from the 2 methods used above
#plot_df <- data.frame(pct = pct, # put PC values into dataframe for plotting
#                      cumu = cumu, 
#                      rank = 1:length(pct))
#ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + # visualize PCs to use in elbow plot
#  geom_text() + 
#  geom_vline(xintercept = 90, color = "grey") + 
#  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#  theme_bw()
PCdims <- 1:pcs # use the minimum PC from the quantitative method above to set the PCs for subsequent steps
#length(PCdims) # how many significant PCs are there?

# Run UMAP dimensionality reduction:
ref <- RunUMAP(ref, dims = PCdims, reduction = "pca") # create UMAP
DimPlot(ref, 
        group.by = "sample.name") # plot by sample ID
DimPlot(ref, 
        group.by = "Region") # plot by region (should only have SmallInt samples)
DimPlot(ref, 
        group.by = "Diagnosis") # plot by diagnosis (should only have healthy adult and pediatric samples)
DimPlot(ref, 
        group.by = "Region.code") # see if cells derive from duodenum, jejunum, or ileum
DimPlot(ref, 
        group.by = "Fraction") # see what type of cell-sorted fractions cells are derived from
DimPlot(ref, 
        group.by = "category") # see broader annotations for all cells; should only see epithelial
DimPlot(ref, 
        group.by = "annotation",
        label = TRUE) # see more specific annotations for all cells

# Save reference as an .h5seurat object:
SaveH5Seurat(ref, filename = "/home/Jayne.Wiarda/scRNAseqEpSI/RefData/HumanSIEpCells_Elmentaite2021_GeneCellFiltered.h5Seurat")
rm(list= ls()[!(ls() %in% c('ref', 'orthoGenes'))]) # clear up space

# Load in the query data:
ep <- LoadH5Seurat('/home/Jayne.Wiarda/scRNAseqEpSI/Seurat/GutEpOnlySubset.h5Seurat')
DefaultAssay(ep) <- 'RNA'

# Filter query data to include only one-to-one gene orthologs & convert to human gene symbols:
querygenes <- as.data.frame(rownames(ep[['RNA']]@counts)) # extract pig gene names from dataset
colnames(querygenes) <- 'gene'
pigGenes <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalAnnot <-gsub("_", "-", pigGenes$FinalAnnot) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[pigGenes$FinalAnnot %in% querygenes$gene, ] # slim down to only genes in our dataset
queryorthos <- intersect(pigGenes$ENSID, orthoGenes$Gene.stable.ID) # find which genes from reference are also one-to-one orthologs
#length(queryorthos) # how many genes are orthologs?
pigGenes <- pigGenes[pigGenes$ENSID %in% queryorthos, ]
#dim(pigGenes)
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$ENSID, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  # retain only unique combinations of pig & human Ensembl IDs, ignoring transcript IDs
#dim(orthoGenes) # should not have same number of rows as in pigGenes
querycounts <- ep[['RNA']]@counts[rownames(ep[['RNA']]@counts) %in% pigGenes$FinalAnnot,]
pigGenes <- pigGenes %>% arrange(factor(FinalAnnot, levels = rownames(querycounts))) # arrange pigGenes in same order as querycounts
orthoGenes <- orthoGenes %>% arrange(factor(Gene.stable.ID, levels = pigGenes$ENSID)) # arrange orthoGenes in same order as pigGenes (and consequently querycounts)
rownames(querycounts) <- orthoGenes$Human.gene.name # change pig genes to human gene names
queryMeta <- ep@meta.data

# Create new Seurat object
query <- CreateSeuratObject(counts = querycounts, # create new Seurat object for our query dataset
                            meta.data = queryMeta)
rm(list= ls()[!(ls() %in% c('ref', 'query', 'ep'))]) # clear up space

# Normalization and integration:
query <- SplitObject(query, split.by = "orig.ident") # split into the original samples that were processed for scRNA-seq
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
  query[[i]] <- RunPCA(query[[i]], # calculate 50 PCs
                       npcs = 50, 
                       verbose = TRUE)
}

# Now save the query data:
saveRDS(query, '/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/HumanComp_EpOnly_query.rds')

## Perform label transfer and mapping:
MappingScores <- list()
CellTypePredictions <- list()
for(i in 1:length(query)) {
  anchors <- FindTransferAnchors(
    reference = ref,
    query = query[[i]],
    reduction = "cca", # opted to use cca since the method is recommended for cross-species mapping 
    dims = 1:30, 
    normalization.method = "LogNormalize",
    k.filter = 132) # reduce k.filter from 200 to 132 since our smallest sample only has 132 cells
  predictions <- TransferData(anchorset = anchors, 
                              refdata = list(cell_type = ref$annotation), # predict query dataset IDs at level of reference data's cluster, lineage, and cell type classifications
                              dims = 1:30,
                              weight.reduction = "cca")
  MapScores <- MappingScore(
    anchors = anchors@anchors,
    combined.object = anchors@object.list[[1]],
    query.neighbors =  slot(object = query[[i]], name = "neighbors")[["query_ref.nn"]],
    query.weights = Tool(object = query[[i]], slot = "TransferData")$weights.matrix,
    query.embeddings = Embeddings(object = query[[i]]),
    ref.embeddings = Embeddings(object = ref),
    nn.method = "annoy",
    # n.trees = n.trees
  )
  MappingScores[[i]] <- MapScores
  CellTypePredictions[[i]] <- predictions
} 

MappingScores <- Reduce(c,MappingScores)
MappingScores <- as.data.frame(MappingScores)
CellTypePredictions <- do.call(rbind, CellTypePredictions)
CellTypePredictions <- as.data.frame(CellTypePredictions)

# Save the mapping & prediction results:
MappingScores$CellBarcodes <- rownames(MappingScores)
CellTypePredictions$CellBarcodes <- rownames(CellTypePredictions)
write_xlsx(MappingScores, '/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/HsSI_EpOnly_Elmentaite2021_MappingScores.xlsx')
write_xlsx(CellTypePredictions, '/home/Jayne.Wiarda/scRNAseqEpSI/MappingPrediction/HsSI_EpOnly_Elmentaite2021_CellTypePredictions.xlsx')

# Incorporate cell prediction & mapping scores into original Seurat object of query data:
ep <- AddMetaData(object = ep, 
                   metadata = c(MappingScores, CellTypePredictions))

FeaturePlot(ep, features = 'MappingScores', 
            reduction = 'tsne', 
            cols = c('yellow', 'red'),
            min.cutoff = 0.6, 
            max.cutoff = 1)

predict <- colnames(ep@meta.data %>% select(starts_with("prediction.score."))) # extract 
names <- sub(".*prediction.score. ", "", predict)   
FeaturePlot(ep,
            features = c(predict),
            reduction = 'tsne',
            ncol = 5) & 
  scale_color_gradientn( colours = c('yellow', 'red'),  limits = c(0, 1), oob = squish) & 
  NoAxes() & NoLegend() 
names

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
#[1] scales_1.1.1          writexl_1.4.0         readxl_1.3.1          dplyr_1.0.7           SeuratDisk_0.0.0.9019 SeuratObject_4.0.2    Seurat_4.0.4          ggplot2_3.3.5        

#loaded via a namespace (and not attached):
#[1] plyr_1.8.6                  igraph_1.2.6                lazyeval_0.2.2              splines_4.1.1               BiocParallel_1.26.2         listenv_0.8.0               scattermore_0.7             GenomeInfoDb_1.28.4         digest_0.6.27              
#[10] htmltools_0.5.2             viridis_0.6.1               fansi_0.5.0                 magrittr_2.0.1              ScaledMatrix_1.0.0          tensor_1.5                  cluster_2.1.2               ROCR_1.0-11                 limma_3.48.3               
#[19] globals_0.14.0              graphlayouts_0.7.1          matrixStats_0.60.1          spatstat.sparse_2.0-0       colorspace_2.0-2            ggrepel_0.9.1               xfun_0.26                   crayon_1.4.1                RCurl_1.98-1.4             
#[28] jsonlite_1.7.2              spatstat.data_2.1-0         survival_3.2-13             zoo_1.8-9                   glue_1.4.2                  polyclip_1.10-0             gtable_0.3.0                zlibbioc_1.38.0             XVector_0.32.0             
#[37] leiden_0.3.9                DelayedArray_0.18.0         BiocSingular_1.8.1          future.apply_1.8.1          SingleCellExperiment_1.14.1 BiocGenerics_0.38.0         abind_1.4-5                 edgeR_3.34.1                DBI_1.1.1                  
#[46] miniUI_0.1.1.1              Rcpp_1.0.7                  viridisLite_0.4.0           xtable_1.8-4                reticulate_1.21             spatstat.core_2.3-0         rsvd_1.0.5                  bit_4.0.4                   htmlwidgets_1.5.4          
#[55] httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.2              ica_1.0-2                   pkgconfig_2.0.3             farver_2.1.0                uwot_0.1.10                 dbplyr_2.1.1                deldir_0.2-10              
#[64] locfit_1.5-9.4              utf8_1.2.2                  labeling_0.4.2              tidyselect_1.1.1            rlang_0.4.11                reshape2_1.4.4              later_1.3.0                 munsell_0.5.0               cellranger_1.1.0           
#[73] tools_4.1.1                 cli_3.0.1                   generics_0.1.0              ggridges_0.5.3              evaluate_0.14               stringr_1.4.0               fastmap_1.1.0               yaml_2.2.1                  goftest_1.2-2              
#[82] knitr_1.34                  bit64_4.0.5                 fitdistrplus_1.1-5          tidygraph_1.2.0             purrr_0.3.4                 RANN_2.6.1                  ggraph_2.0.5                pbapply_1.5-0               future_1.22.1              
#[91] nlme_3.1-152                mime_0.11                   hdf5r_1.3.4                 compiler_4.1.1              beeswarm_0.4.0              plotly_4.9.4.1              png_0.1-7                   spatstat.utils_2.2-0        tibble_3.1.4               
#[100] tweenr_1.0.2                stringi_1.7.4               lattice_0.20-44             Matrix_1.3-4                vctrs_0.3.8                 pillar_1.6.2                lifecycle_1.0.0             spatstat.geom_2.2-2         lmtest_0.9-38              
#[109] RcppAnnoy_0.0.19            BiocNeighbors_1.10.0        data.table_1.14.0           cowplot_1.1.1               bitops_1.0-7                irlba_2.3.3                 httpuv_1.6.3                patchwork_1.1.1             GenomicRanges_1.44.0       
#[118] R6_2.5.1                    promises_1.2.0.1            KernSmooth_2.23-20          gridExtra_2.3               vipor_0.4.5                 IRanges_2.26.0              parallelly_1.28.1           codetools_0.2-18            gtools_3.9.2               
#[127] MASS_7.3-54                 assertthat_0.2.1            SummarizedExperiment_1.22.0 withr_2.4.2                 sctransform_0.3.2           S4Vectors_0.30.0            GenomeInfoDbData_1.2.6      mgcv_1.8-36                 miloR_1.0.0                
#[136] beachmat_2.8.1              grid_4.1.1                  rpart_4.1-15                tidyr_1.1.3                 rmarkdown_2.11              MatrixGenerics_1.4.3        Rtsne_0.15                  ggforce_0.3.3               Biobase_2.52.0             
#[145] shiny_1.6.0                 ggbeeswarm_0.6.0       