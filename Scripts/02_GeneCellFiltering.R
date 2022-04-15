# Title: 02 - Gene and Cell Filtering
# Date: 2021Sept20
# Author: Jayne Wiarda

library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(matrixStats) 
library(writexl)
library(DropletUtils)
library(readxl)

## Create Seurat object:
data_dir <- c(DUOD2 = "/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/DUOD2strainedCounts", # specify paths to output files from each sample that have been strained with SoupX; listed alphabetically
              JEJ2 = "/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/JEJ2strainedCounts",               
              IPP2 = "/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/IPP2strainedCounts", 
              NoPP2 = "/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/NoPP2strainedCounts")
library_id <- c("DUOD2", "JEJ2", "IPP2", "NoPP2") # set the Sample IDs to use (in same order as in data_dir)
lapply(data_dir, dir) # Should show barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx.gz for each sample listed
scRNA_data <- Read10X(data.dir = data_dir) # read the 10X data from all samples into a data matrix
seurat_object = CreateSeuratObject(counts = scRNA_data) # create a Seurat object of the data matrix
#seurat_object # see how many genes/cells present in Seurat object

## Create a dataframe with cell phenotype data information (pDat):
pDat <-data.frame(barcode = colnames(seurat_object)) # create a dataframe of the cell barcodes
pDat$SampleID <- gsub("([^_]+).+", "\\1", pDat$barcode, perl = TRUE) # make column of SampleIDs for each row entry/cell

## Make a column of barcodes corresponding to .cloupe file from Cell Ranger outputs:
pDat$BarBak <- pDat$barcode # make new column called BarBak with same entries as barcode column
pDat <- pDat %>% separate(BarBak, c("Sam","Loupe")) # separate BarBak entries into sample ID ("Sam") and barcode ID ("Loupe")
pDat <- pDat[,-3] # remove Sam column
for (i in seq_along(library_id)){ # Add -# suffix to end of Loupe column entries, corresponding to sample order in library_id. We will need these to correspond with barcodes in the .cloupe file from CellRanger outputs in order to later edit data for Loupe Cell Browser
  pDat$Loupe <- ifelse(pDat$SampleID == library_id[i], paste0(pDat$Loupe, paste0("-",i)), pDat$Loupe)
}
rownames(pDat) <- pDat$Loupe # make Loupe barcodes the rownames of pDat

## Check the pDat format:
#tail(pDat, n = 3) # pDat should look something like this

table(duplicated(pDat$Loupe)) # make sure none of our Loupe barcode IDs are duplicated

## Check the number of cells per sample:
table(pDat$SampleID)

## Create count matrix (cDat)
cDat <-  as.matrix(GetAssayData(object = seurat_object, slot = 'counts')) # create a new counts matrix; columns = cells & rows = genes
#dim(cDat)	# obtain matrix dimensions indicating number of genes x number of cells:

## Filter out lowly expressed genes from the count matrix
keep <- rowSums(cDat) > 0 # specify rows with genes that have no expression across all cells of the dataset
cDat <- cDat[keep,]  # keep only genes that have > 0 transcript total across all cells of the dataset
#dim (cDat)	# obtain new matrix dimensions after gene filtering indicating number of genes x number of cells:

## Create a dataframe with feature data (fDat):
fDat <- data.frame(ID = rownames(cDat)) # create a dataframe of the filtered genes
rownames(fDat) <- fDat$ID # make genes the row names of dataframe
#head(fDat, n = 3) # dataframe should look like this:

## Add gene symbol information to fDat:
n_last <- 18 # Ensembl IDs are 18 characters long
fDat$EnsemblID <- substr(fDat$ID, nchar(fDat$ID) - n_last + 1, nchar(fDat$ID)) # create a column containing only Ensembl IDs in fDat

## Bring in updated gene symbol annotations:
geneannot <- read_excel('/home/Jayne.Wiarda/scRNAseqEpSI/UpdatedGeneNameListForSus97GTF.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
geneannot <- filter(geneannot, ENSID %in% fDat$EnsemblID)# reduce new gene sybmol annotation to include only genes in the filtered gene list for our single-cell dataset
dim(geneannot)

## Incorporate updated gene symbol annotations into fDat:
fDat <- merge(fDat, geneannot, by.x = "EnsemblID", by.y = "ENSID")
#head(fDat)

## Extract mitochondrial gene information and add to fDat:
# Requires download and unzipping of Sus_scrofa.Sscrofa11.1.97.gtf.gz annotation
gtf <- "//home/Jayne.Wiarda/scRNAseqIleumAtlas/SS_annotation/Sus_scrofa.Sscrofa11.1.97.gtf" # specify file path to Sus scrofa 11.1 version 100 annotation file
mitoGenes <-  system2("grep", args = c('^MT', gtf, "| grep -o 'ENSSSCG[0-9]*' | uniq"), stdout = TRUE) # extract mitochondrial gene information from annotation file
fDat$Mitochondrial <- fDat$EnsemblID %in% mitoGenes               # add column of mitochondrial gene information to feature dataframe
length(mitoGenes) # see that we have 37 mitochondrial genes:
# [1] 37
#table(fDat$Mitochondrial) # see that we have the same amount of "TRUE" mito genes as in our mitoGenes list (37 genes)

## Modify count matrix (cDat) barcodes & gene names:
fDat <- fDat %>% arrange(factor(ID, levels = rownames(cDat))) # arrange fDat genes in same order as in cDat
all(rownames(cDat) == fDat$ID) # make sure genes that are rownames in cDat are in same exact order in fDat rows
rownames(cDat) <- fDat$FinalAnnot # change cDat rownames to fDat updated gene symbol annotations
#head(rownames(cDat), n = 3) # see that changes occurred
length(unique(rownames(cDat))) # make sure the number of unique row names in cDat is equal to the total number of rows in cDat
nrow(cDat)

########################################
#### Filter out poor quality cells: ####
########################################

## Assess single-cell sequencing depths and number of genes detected:
pDat$UmiSums<- colSums(cDat) # find the total unique reads detected per cell
pDat$GenesDetected <- colSums(cDat!=0) # find the total number of genes detected per cell
ggplot(pDat, aes(x=SampleID, y=GenesDetected, fill= SampleID)) + # create violin plot of number of total genes detected per cell in each sample
  geom_violin(draw_quantiles=0.5)+
  ylab("Total number of genes detected")
ggplot(pDat, aes(x=SampleID,y=UmiSums, fill=SampleID)) + # create violin plot of number of unique reads (UMIs) detected per cell in each sample
  geom_violin(draw_quantiles=0.5)+
  ylab("Total number of molecules(depth)")

## Plot 50 most highly expressed genes:
mRel <- t(t(cDat)/colSums(cDat)) # calculate percentage of total expression for each gene in each cell
table(colSums(mRel)) # should all total 1
rownames(mRel)  <- fDat$ID # assign rownames
topExpressed <- rowMedians(mRel) # find median expression level to see library size normalized expression
names(topExpressed) <- rownames(mRel) # assign corresponding gene names to median expression values
topExpressed <- topExpressed %>% sort(.,decreasing=TRUE) %>% names # sort gene names according to highest to lowest median library-normalized expression
plotData <- t(mRel)[,topExpressed[1:50]] %>% reshape2::melt() %>% # create dataframe of top 50 expressed genes individual cell expression levels
  dplyr::rename(Cell=Var1, Gene=Var2, RelativeExpression=value)
plotData <- merge(plotData, fDat, by.x = "Gene", by.y = "ID", all.x = TRUE) # Add mitochondrial gene information for top 50 expressed genes
ggplot(plotData, aes(x=Gene, y=RelativeExpression, color= Mitochondrial)) + # create plot of top 50 expressed genes and whether they are mitochondrial genes... should see lots of ribosomal and some mitochondrial genes
  geom_boxplot() +     
  coord_flip() +     
  theme_bw()

## Plot 50 most frequently expressed genes:
freqOfExp <- cDat!=0 # list only genes with non-zero expression and how many cells express the gene
rownames(freqOfExp) <- fDat$ID # make a table of whether or not a gene is expressed in each cell
freqOfExp <- sort(rowSums(freqOfExp)/ncol(freqOfExp),decreasing=TRUE) # list percentages of cells expressing a gene, from highest to lowest
plotData <- data.frame("Gene"=names(freqOfExp),"Frequency"=freqOfExp) # create dataframe of gene names and frequencies
plotData <- merge(plotData, fDat, by.x = "Gene", by.y = "ID", all.x = TRUE, sorted =FALSE) # add mitochodrial gene information
plotData <- plotData[order(plotData$Frequency, decreasing= TRUE), ] # sort from highest to lowest frequency

## most frequently detected genes across all cells
ggplot(plotData[1:50,], aes(x=factor(Gene,levels=Gene), y=Frequency, color= Mitochondrial)) + # plot 50 most frequently expressed genes... should see lots of ribosomal and some mitochondrial genes
  geom_bar(stat="identity", fill ="white") +
  coord_flip() +     
  xlab("Gene") +     
  theme_bw()

## Assess cell quality:
theme_set(theme_grey()) # set plot theme
mMito <- cDat[fDat$Mitochondrial,] # extract read count of only mitochondrial genes
idtop <- fDat[fDat$Name %in% names(freqOfExp)[1:50],"ID"] # get IDs of top 50 most frequently expressed genes
mTop <- cDat[idtop,]!=0 # Print TRUE/FALSE expression of top 50 genes alone. rows - top 50 genes expression, column - all cells.
pDat$prcntTop <- colSums(mTop)/50 # calculate percentage/frequency of expression of top 50 genes in all cells. colsum is sum of reads on these top50 genes in each cells (column wise) devided by 50 (50)
pDat$prcntMito <- colSums(mMito)/colSums(cDat) # add column of % mitochondrial gene expression for each cell
ggplot(pDat, aes(x=prcntMito, y=GenesDetected, color = SampleID))+ # plot % mito genes vs # genes detected per cell
  geom_point() + 
  facet_wrap(~SampleID, nrow =1)+
  theme_get() + 
  ylab("#Genes detected per cell")

## Check for barcode duplications:
pDat$DuplicatedBarcodes <- duplicated(rownames(pDat)) | duplicated(rownames(pDat), fromLast = TRUE)
table(pDat$DuplicatedBarcodes) # since no barcodes are duplicated, we won't have to remove any cells due to this concern, but we will still leave this criteria in our cell filtering parameters

## Establish QC thresholds:
# Look at histograms to determine if there are obvious cutoff values:
ggplot(pDat, aes(x=prcntMito,y=..density..)) + # plot % mitochondrial reads to find good cutoff value for QC filtering
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, .5, .025), lim = c(0, .5)) + 
  ylim(0,30) +
  facet_wrap(~SampleID) +
  geom_vline(aes(xintercept=.2),color="red",lty="longdash") + # move this cutoff line where you see fit; 20% mitochondrial reads seems like a good cutoff
  RotatedAxis()
ggplot(pDat, aes(x=GenesDetected,y=..density..)) + # plot # genes detected per cell to find good cutoff value for QC filtering
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 2000, 250), lim = c(0, 2000)) + 
  RotatedAxis() +
  geom_vline(aes(xintercept=500),color="red",lty="longdash") + # move this cutoff line where you see fit; 500 genes detected seems like a good cutoff
  facet_wrap(~SampleID) 
ggplot(pDat, aes(x=UmiSums,y=..density..)) + # plot # UMIs detected per cell to find good cutoff value for QC filtering
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 5000, 250), lim = c(0, 5000)) + 
  RotatedAxis() +
  geom_vline(aes(xintercept=1250),color="red",lty="longdash") + # move this cutoff line where you see fit; 1250 UMIs seems like a good cutoff
  facet_wrap(~SampleID)

pDat <- mutate(pDat, 
               PassViability=prcntMito < 0.2, # only include cells with total mitochondrial reads under 12.5%
               PassGenesDet=GenesDetected > 500, # only consider cells with total genes detected more than 550
               PassLibSize=UmiSums > 1250, # only consider cells with greater than 1250 total UMIs
               PassBarcodeFreq=DuplicatedBarcodes==FALSE, # only consider cells with non-repeated barcodes
               PassAll= PassViability & PassGenesDet & PassLibSize & PassBarcodeFreq) # list whether or not cells pass all filtering criteria
rownames(pDat) <- pDat$Loupe # make sure the pDat rownames correspond to Loupe barcodes!

## Overview over cells removed
table(pDat$SampleID,pDat$PassGenesDet) # how many cells had at least 500 genes detected?
table(pDat$SampleID,pDat$PassLibSize) # how many cells had at least 1250 UMIs detected?
table(pDat$SampleID,pDat$PassViability)  # how many cells had less than 20% mitochondrial reads?
table(pDat$SampleID,pDat$PassBarcodeFreq)  # how many cells had non-repeated barcodes?
table(pDat$SampleID,pDat$PassAll) # how many cells passed all filtering criteria?

## Save Data
stopifnot(identical(as.character(fDat$FinalAnnot),rownames(cDat)))
out <- list()
out[["counts"]] <- cDat
out[["phenoData"]] <- pDat
out[["featureData"]] <- fDat
saveRDS(out,file=file.path("/home/Jayne.Wiarda/scRNAseqEpSI/QC/", "GutEpQC.rds")) # this saves all of our information before filtering out low quality cells
write.table(fDat, 
            file = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/GeneInfo.txt") # export the feature information separately...this will come in handy later
write_xlsx(x = fDat, 
           path = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/GeneInfo.xlsx", # export the feature information separately...this will come in handy later
           col_names = TRUE) # this is the file used for cluster-specific cell type 
write.table(ssc_genes, 
            file = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/UnfilteredGeneInfo.txt") # export the feature information separately...this will come in handy later
write_xlsx(x = ssc_genes, 
           path = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/UnfilteredGeneInfo.xlsx", # export the feature information separately...this will come in handy later
           col_names = TRUE) # this is the file used for cluster-specific cell type 

## Filter out poor quality cells:
dim(cDat) 
cDat <- cDat[,pDat$PassAll] # keep only cells that passed all criteria in QC
dim(cDat) 
pDat <- pDat[pDat$PassAll,]
dim(pDat) # make sure rows here matches columns in cDat

## subset data & save filtered data in CellRanger formats:
All <- CreateSeuratObject(counts = cDat) # create Seurat object of counts & pheno data
All$SampleID <- substr(colnames(cDat), 1, 3) 
Idents(All) <- "SampleID"

DUOD2 <- subset(All, ident = "DUO") # subset only cells originating from DUOD2 sample
write10xCounts(x = DUOD2@assays$RNA@counts, path = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/DUOD2onlyFilteredQC", version = "3") # create CellRanger-like output files of our subsetted Seurat object
rm(DUOD2)

JEJ2 <- subset(All, ident = "JEJ")
write10xCounts(x = JEJ2@assays$RNA@counts, path = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/JEJ2onlyFilteredQC", version = "3")
rm(JEJ2)

IPP2 <- subset(All, ident = "IPP")
write10xCounts(x = IPP2@assays$RNA@counts, path = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/IPP2onlyFilteredQC", version = "3")
rm(IPP2)

NoPP2 <- subset(All, ident = "NoP")
write10xCounts(x = NoPP2@assays$RNA@counts, path = "/home/Jayne.Wiarda/scRNAseqEpSI/QC/NoPP2onlyFilteredQC", version = "3")
rm(NoPP2)

rm(All)

## In total, we have generated count matrices for each individual sample that have strained (SoupX) and filtered (genes with no expression removed) gene lists, as well as filtering to remove poor quality cells
## At this point, we transfer our _SampleNameHere_onlyFilteredQC directories and accompanying files over to /project/fsepru/PIGscRNAseq1/100/ and run first part of doublet removal code using Python script

#sessionInfo()
#R version 4.1.1 (2021-08-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] readxl_1.3.1                DropletUtils_1.12.2         SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0 Biobase_2.52.0              GenomicRanges_1.44.0        GenomeInfoDb_1.28.4         IRanges_2.26.0              S4Vectors_0.30.0            BiocGenerics_0.38.0        
#[11] MatrixGenerics_1.4.3        writexl_1.4.0               matrixStats_0.60.1          tidyr_1.1.3                 dplyr_1.0.7                 SeuratObject_4.0.2          Seurat_4.0.4                ggplot2_3.3.5              

#loaded via a namespace (and not attached):
#[1] plyr_1.8.6                igraph_1.2.6              lazyeval_0.2.2            splines_4.1.1             BiocParallel_1.26.2       listenv_0.8.0             scattermore_0.7           digest_0.6.27             htmltools_0.5.2           viridis_0.6.1             fansi_0.5.0              
#[12] magrittr_2.0.1            ScaledMatrix_1.0.0        tensor_1.5                cluster_2.1.2             ROCR_1.0-11               limma_3.48.3              graphlayouts_0.7.1        globals_0.14.0            R.utils_2.10.1            spatstat.sparse_2.0-0     colorspace_2.0-2         
#[23] ggrepel_0.9.1             xfun_0.26                 crayon_1.4.1              RCurl_1.98-1.4            jsonlite_1.7.2            spatstat.data_2.1-0       survival_3.2-13           zoo_1.8-9                 glue_1.4.2                polyclip_1.10-0           gtable_0.3.0             
#[34] zlibbioc_1.38.0           XVector_0.32.0            leiden_0.3.9              DelayedArray_0.18.0       BiocSingular_1.8.1        Rhdf5lib_1.14.2           future.apply_1.8.1        HDF5Array_1.20.0          abind_1.4-5               scales_1.1.1              DBI_1.1.1                
#[45] edgeR_3.34.1              miniUI_0.1.1.1            Rcpp_1.0.7                viridisLite_0.4.0         xtable_1.8-4              reticulate_1.21           spatstat.core_2.3-0       dqrng_0.3.0               rsvd_1.0.5                htmlwidgets_1.5.4         httr_1.4.2               
#[56] RColorBrewer_1.1-2        ellipsis_0.3.2            ica_1.0-2                 farver_2.1.0              scuttle_1.2.1             pkgconfig_2.0.3           R.methodsS3_1.8.1         uwot_0.1.10               dbplyr_2.1.1              deldir_0.2-10             locfit_1.5-9.4           
#[67] utf8_1.2.2                labeling_0.4.2            tidyselect_1.1.1          rlang_0.4.11              reshape2_1.4.4            later_1.3.0               cellranger_1.1.0          munsell_0.5.0             tools_4.1.1               cli_3.0.1                 generics_0.1.0           
#[78] ggridges_0.5.3            evaluate_0.14             stringr_1.4.0             fastmap_1.1.0             yaml_2.2.1                goftest_1.2-2             knitr_1.34                tidygraph_1.2.0           fitdistrplus_1.1-5        purrr_0.3.4               RANN_2.6.1               
#[89] ggraph_2.0.5              pbapply_1.5-0             future_1.22.1             nlme_3.1-152              sparseMatrixStats_1.4.2   mime_0.11                 R.oo_1.24.0               rstudioapi_0.13           compiler_4.1.1            beeswarm_0.4.0            plotly_4.9.4.1           
#[100] png_0.1-7                 spatstat.utils_2.2-0      tweenr_1.0.2              tibble_3.1.4              stringi_1.7.4             lattice_0.20-44           Matrix_1.3-4              permute_0.9-5             vegan_2.5-7               vctrs_0.3.8               pillar_1.6.2             
#[111] lifecycle_1.0.0           rhdf5filters_1.4.0        spatstat.geom_2.2-2       lmtest_0.9-38             BiocNeighbors_1.10.0      RcppAnnoy_0.0.19          data.table_1.14.0         cowplot_1.1.1             bitops_1.0-7              irlba_2.3.3               httpuv_1.6.3             
#[122] patchwork_1.1.1           R6_2.5.1                  promises_1.2.0.1          KernSmooth_2.23-20        gridExtra_2.3             vipor_0.4.5               parallelly_1.28.1         codetools_0.2-18          gtools_3.9.2              MASS_7.3-54               assertthat_0.2.1         
#[133] rhdf5_2.36.0              withr_2.4.2               sctransform_0.3.2         GenomeInfoDbData_1.2.6    miloR_1.0.0               mgcv_1.8-36               grid_4.1.1                rpart_4.1-15              beachmat_2.8.1            rmarkdown_2.11            DelayedMatrixStats_1.14.3
#[144] Rtsne_0.15                ggforce_0.3.3             shiny_1.6.0               ggbeeswarm_0.6.0