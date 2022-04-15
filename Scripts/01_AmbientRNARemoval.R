# Title: 01 - Ambient RNA Removal
# Date: 2021Sept20
# Author: Jayne Wiarda

library(SoupX)
library(DropletUtils)
library(ggplot2)

sc = load10X('/home/Jayne.Wiarda/scRNAseqEpSI/CellRanger/Pigdataset2/DUOD2') # Provide file path to folder containing folders of both raw and filtered Cell Ranger outputs to create a SoupChannel object
#head(sc$metaData, n = 3) # double check metadata for cluster ID and tsne was read in

# Plot the data:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot; this should be identical to what would be brought up in Loupe Cell Browser using .cloupe file from Cell Ranger outputs for just that single sample (not the same coordinates for that sample in the .cloupe file that contains all of the samples)
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 

# Now let's look at some specific gene expression:
#dd$CD3E = sc$toc["CD3E_ENSSSCG00000040140", ] # make column of gene expression values for CD3E
#dd$IgLambdaV = sc$toc["ENSSSCG00000038719", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$CD79B = sc$toc["CD79B_ENSSSCG00000017283", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$FABP6 = sc$toc["FABP6_ENSSSCG00000017037", ] # make column of gene expression values for gene that codes for FABP6
#dd$EPCAM = sc$toc["EPCAM_ENSSSCG00000008429", ] # make column of gene expression values for gene that codes for EPCAM
#dd$GNLY = sc$toc["GNLY_ENSSSCG00000008228", ] # make column of gene expression values for gene that codes for GNLY
#dd$HBB = sc$toc["HBB_ENSSSCG00000014725", ] # make column of gene expression values for gene that codes for GNLY
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) # which cells express this gene?
#plotMarkerMap(sc, "CD3E_ENSSSCG00000040140") # if we assumed all cells were nothing but soup, which cells still show higher than expected expression for the gene (TRUE = expression levels higher than expected if cell was just soup, so likely real expression)? This just gives us an idea of soup expression, this is NOT a formal analysis used for removing the soup RNA.
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = IgLambdaV > 0))
#plotMarkerMap(sc, "ENSSSCG00000038719")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0))
#plotMarkerMap(sc, "CD79B_ENSSSCG00000017283")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FABP6 > 0)) 
#plotMarkerMap(sc, "FABP6_ENSSSCG00000017037")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = EPCAM > 0)) 
#plotMarkerMap(sc, "EPCAM_ENSSSCG00000008429")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = GNLY > 0))
#plotMarkerMap(sc, "GNLY_ENSSSCG00000008228")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0))
#plotMarkerMap(sc, "HBB_ENSSSCG00000014725")
# What we see is some misplaced gene expression, indicating we have RNA soup and need to remove....

## Estimate the contaminated RNA fraction automatically:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

## Now let's see which genes were most affected by our correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
#tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
#tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

## Let's make sure the results for some of our genes make sense by visualizing on a plot:
#plotChangeMap(sc, out, "CD3E_ENSSSCG00000040140") # see which cells had gene expression mainly related to soup for various genes
#plotChangeMap(sc, out, "ENSSSCG00000038719")
#plotChangeMap(sc, out, "CD79B_ENSSSCG00000017283")
#plotChangeMap(sc, out, "FABP6_ENSSSCG00000017037")
#plotChangeMap(sc, out, "EPCAM_ENSSSCG00000008429")
#plotChangeMap(sc, out, "GNLY_ENSSSCG00000008228")
#plotChangeMap(sc, out, "HBB_ENSSSCG00000014725")

## We are happy with our results and need to save our strained count matrices in a new location:
write10xCounts("/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/DUOD2strainedCounts", out, version = "3")
rm(dd, mids, out, sc, cntSoggy, cntStrained)

## Now let's perform similar analyses on the rest of our samples:

#JEJ2:
sc = load10X('/home/Jayne.Wiarda/scRNAseqEpSI/CellRanger/Pigdataset2/JEJ2') # Provide file path to folder containing folders of both raw and filtered Cell Ranger outputs to create a SoupChannel object
#head(sc$metaData, n = 3) # double check metadata for cluster ID and tsne was read in

# Plot the data:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot; this should be identical to what would be brought up in Loupe Cell Browser using .cloupe file from Cell Ranger outputs for just that single sample (not the same coordinates for that sample in the .cloupe file that contains all of the samples)
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 

# Now let's look at some specific gene expression:
#dd$CD3E = sc$toc["CD3E_ENSSSCG00000040140", ] # make column of gene expression values for CD3E
#dd$IgLambdaV = sc$toc["ENSSSCG00000038719", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$CD79B = sc$toc["CD79B_ENSSSCG00000017283", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$FABP6 = sc$toc["FABP6_ENSSSCG00000017037", ] # make column of gene expression values for gene that codes for FABP6
#dd$EPCAM = sc$toc["EPCAM_ENSSSCG00000008429", ] # make column of gene expression values for gene that codes for EPCAM
#dd$GNLY = sc$toc["GNLY_ENSSSCG00000008228", ] # make column of gene expression values for gene that codes for GNLY
#dd$HBB = sc$toc["HBB_ENSSSCG00000014725", ] # make column of gene expression values for gene that codes for GNLY
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) # which cells express this gene?
#plotMarkerMap(sc, "CD3E_ENSSSCG00000040140") # if we assumed all cells were nothing but soup, which cells still show higher than expected expression for the gene (TRUE = expression levels higher than expected if cell was just soup, so likely real expression)? This just gives us an idea of soup expression, this is NOT a formal analysis used for removing the soup RNA.
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = IgLambdaV > 0))
#plotMarkerMap(sc, "ENSSSCG00000038719")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0))
#plotMarkerMap(sc, "CD79B_ENSSSCG00000017283")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FABP6 > 0)) 
#plotMarkerMap(sc, "FABP6_ENSSSCG00000017037")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = EPCAM > 0)) 
#plotMarkerMap(sc, "EPCAM_ENSSSCG00000008429")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = GNLY > 0))
#plotMarkerMap(sc, "GNLY_ENSSSCG00000008228")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0))
#plotMarkerMap(sc, "HBB_ENSSSCG00000014725")
# What we see is some misplaced gene expression, indicating we have RNA soup and need to remove....

## Estimate the contaminated RNA fraction automatically:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

## Now let's see which genes were most affected by our correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
#tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
#tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

## Let's make sure the results for some of our genes make sense by visualizing on a plot:
#plotChangeMap(sc, out, "CD3E_ENSSSCG00000040140") # see which cells had gene expression mainly related to soup for various genes
#plotChangeMap(sc, out, "ENSSSCG00000038719")
#plotChangeMap(sc, out, "CD79B_ENSSSCG00000017283")
#plotChangeMap(sc, out, "FABP6_ENSSSCG00000017037")
#plotChangeMap(sc, out, "EPCAM_ENSSSCG00000008429")
#plotChangeMap(sc, out, "GNLY_ENSSSCG00000008228")
#plotChangeMap(sc, out, "HBB_ENSSSCG00000014725")

## We are happy with our results and need to save our strained count matrices in a new location:
write10xCounts("/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/JEJ2strainedCounts", out, version = "3")
rm(dd, mids, out, sc, cntSoggy, cntStrained)

#IPP2:
sc = load10X('/home/Jayne.Wiarda/scRNAseqEpSI/CellRanger/Pigdataset2/IPP2') # Provide file path to folder containing folders of both raw and filtered Cell Ranger outputs to create a SoupChannel object
#head(sc$metaData, n = 3) # double check metadata for cluster ID and tsne was read in

# Plot the data:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot; this should be identical to what would be brought up in Loupe Cell Browser using .cloupe file from Cell Ranger outputs for just that single sample (not the same coordinates for that sample in the .cloupe file that contains all of the samples)
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 

# Now let's look at some specific gene expression:
#dd$CD3E = sc$toc["CD3E_ENSSSCG00000040140", ] # make column of gene expression values for CD3E
#dd$IgLambdaV = sc$toc["ENSSSCG00000038719", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$CD79B = sc$toc["CD79B_ENSSSCG00000017283", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$FABP6 = sc$toc["FABP6_ENSSSCG00000017037", ] # make column of gene expression values for gene that codes for FABP6
#dd$EPCAM = sc$toc["EPCAM_ENSSSCG00000008429", ] # make column of gene expression values for gene that codes for EPCAM
#dd$GNLY = sc$toc["GNLY_ENSSSCG00000008228", ] # make column of gene expression values for gene that codes for GNLY
#dd$HBB = sc$toc["HBB_ENSSSCG00000014725", ] # make column of gene expression values for gene that codes for GNLY
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) # which cells express this gene?
#plotMarkerMap(sc, "CD3E_ENSSSCG00000040140") # if we assumed all cells were nothing but soup, which cells still show higher than expected expression for the gene (TRUE = expression levels higher than expected if cell was just soup, so likely real expression)? This just gives us an idea of soup expression, this is NOT a formal analysis used for removing the soup RNA.
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = IgLambdaV > 0))
#plotMarkerMap(sc, "ENSSSCG00000038719")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0))
#plotMarkerMap(sc, "CD79B_ENSSSCG00000017283")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FABP6 > 0)) 
#plotMarkerMap(sc, "FABP6_ENSSSCG00000017037")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = EPCAM > 0)) 
#plotMarkerMap(sc, "EPCAM_ENSSSCG00000008429")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = GNLY > 0))
#plotMarkerMap(sc, "GNLY_ENSSSCG00000008228")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0))
#plotMarkerMap(sc, "HBB_ENSSSCG00000014725")
# What we see is some misplaced gene expression, indicating we have RNA soup and need to remove....

## Estimate the contaminated RNA fraction automatically:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

## Now let's see which genes were most affected by our correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
#tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
#tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

## Let's make sure the results for some of our genes make sense by visualizing on a plot:
#plotChangeMap(sc, out, "CD3E_ENSSSCG00000040140") # see which cells had gene expression mainly related to soup for various genes
#plotChangeMap(sc, out, "ENSSSCG00000038719")
#plotChangeMap(sc, out, "CD79B_ENSSSCG00000017283")
#plotChangeMap(sc, out, "FABP6_ENSSSCG00000017037")
#plotChangeMap(sc, out, "EPCAM_ENSSSCG00000008429")
#plotChangeMap(sc, out, "GNLY_ENSSSCG00000008228")
#plotChangeMap(sc, out, "HBB_ENSSSCG00000014725")

## We are happy with our results and need to save our strained count matrices in a new location:
write10xCounts("/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/IPP2strainedCounts", out, version = "3")
rm(dd, mids, out, sc, cntSoggy, cntStrained)

#NoPP2:
sc = load10X('/home/Jayne.Wiarda/scRNAseqEpSI/CellRanger/Pigdataset2/NoPP2') # Provide file path to folder containing folders of both raw and filtered Cell Ranger outputs to create a SoupChannel object
#head(sc$metaData, n = 3) # double check metadata for cluster ID and tsne was read in

# Plot the data:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot; this should be identical to what would be brought up in Loupe Cell Browser using .cloupe file from Cell Ranger outputs for just that single sample (not the same coordinates for that sample in the .cloupe file that contains all of the samples)
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 

# Now let's look at some specific gene expression:
#dd$CD3E = sc$toc["CD3E_ENSSSCG00000040140", ] # make column of gene expression values for CD3E
#dd$IgLambdaV = sc$toc["ENSSSCG00000038719", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$CD79B = sc$toc["CD79B_ENSSSCG00000017283", ] # make column of gene expression values for gene that codes for Ig lambda V region
#dd$FABP6 = sc$toc["FABP6_ENSSSCG00000017037", ] # make column of gene expression values for gene that codes for FABP6
#dd$EPCAM = sc$toc["EPCAM_ENSSSCG00000008429", ] # make column of gene expression values for gene that codes for EPCAM
#dd$GNLY = sc$toc["GNLY_ENSSSCG00000008228", ] # make column of gene expression values for gene that codes for GNLY
#dd$HBB = sc$toc["HBB_ENSSSCG00000014725", ] # make column of gene expression values for gene that codes for GNLY
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) # which cells express this gene?
#plotMarkerMap(sc, "CD3E_ENSSSCG00000040140") # if we assumed all cells were nothing but soup, which cells still show higher than expected expression for the gene (TRUE = expression levels higher than expected if cell was just soup, so likely real expression)? This just gives us an idea of soup expression, this is NOT a formal analysis used for removing the soup RNA.
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = IgLambdaV > 0))
#plotMarkerMap(sc, "ENSSSCG00000038719")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0))
#plotMarkerMap(sc, "CD79B_ENSSSCG00000017283")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FABP6 > 0)) 
#plotMarkerMap(sc, "FABP6_ENSSSCG00000017037")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = EPCAM > 0)) 
#plotMarkerMap(sc, "EPCAM_ENSSSCG00000008429")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = GNLY > 0))
#plotMarkerMap(sc, "GNLY_ENSSSCG00000008228")
#ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0))
#plotMarkerMap(sc, "HBB_ENSSSCG00000014725")
# What we see is some misplaced gene expression, indicating we have RNA soup and need to remove....

## Estimate the contaminated RNA fraction automatically:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

## Now let's see which genes were most affected by our correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
#tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
#tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

## Let's make sure the results for some of our genes make sense by visualizing on a plot:
#plotChangeMap(sc, out, "CD3E_ENSSSCG00000040140") # see which cells had gene expression mainly related to soup for various genes
#plotChangeMap(sc, out, "ENSSSCG00000038719")
#plotChangeMap(sc, out, "CD79B_ENSSSCG00000017283")
#plotChangeMap(sc, out, "FABP6_ENSSSCG00000017037")
#plotChangeMap(sc, out, "EPCAM_ENSSSCG00000008429")
#plotChangeMap(sc, out, "GNLY_ENSSSCG00000008228")
#plotChangeMap(sc, out, "HBB_ENSSSCG00000014725")

## We are happy with our results and need to save our strained count matrices in a new location:
write10xCounts("/home/Jayne.Wiarda/scRNAseqEpSI/SoupX/NoPP2strainedCounts", out, version = "3")
rm(dd, mids, out, sc, cntSoggy, cntStrained)

## In total, we have generated 'strained' count matrices that have removed lots of the ambient RNA (soup) within each sample
## Note that the new counts are not necessarily integers anymore due to the straining method used. If desired (but not used here), counts can still be rounded to the nearest whole number

#sessionInfo()
#R version 4.1.1 (2021-08-10)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS

#Matrix products: default
#BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
#[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8
#[12] LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggplot2_3.3.5               DropletUtils_1.12.2         SingleCellExperiment_1.14.1 SummarizedExperiment_1.22.0 Biobase_2.52.0              GenomicRanges_1.44.0        GenomeInfoDb_1.28.4         IRanges_2.26.0              S4Vectors_0.30.0            BiocGenerics_0.38.0         MatrixGenerics_1.4.3       
#[12] matrixStats_0.60.1          SoupX_1.5.2                

#loaded via a namespace (and not attached):
#[1] utf8_1.2.2                reticulate_1.21           R.utils_2.10.1            tidyselect_1.1.1          htmlwidgets_1.5.4         grid_4.1.1                BiocParallel_1.26.2       Rtsne_0.15                munsell_0.5.0             ScaledMatrix_1.0.0        codetools_0.2-18         
#[12] ica_1.0-2                 future_1.22.1             miniUI_0.1.1.1            withr_2.4.2               colorspace_2.0-2          knitr_1.34                Seurat_4.0.4              ROCR_1.0-11               tensor_1.5                listenv_0.8.0             labeling_0.4.2           
#[23] GenomeInfoDbData_1.2.6    polyclip_1.10-0           farver_2.1.0              rhdf5_2.36.0              parallelly_1.28.1         vctrs_0.3.8               generics_0.1.0            xfun_0.26                 R6_2.5.1                  ggbeeswarm_0.6.0          graphlayouts_0.7.1       
#[34] rsvd_1.0.5                locfit_1.5-9.4            miloR_1.0.0               bitops_1.0-7              rhdf5filters_1.4.0        spatstat.utils_2.2-0      DelayedArray_0.18.0       assertthat_0.2.1          promises_1.2.0.1          scales_1.1.1              ggraph_2.0.5             
#[45] beeswarm_0.4.0            gtable_0.3.0              beachmat_2.8.1            globals_0.14.0            goftest_1.2-2             tidygraph_1.2.0           rlang_0.4.11              splines_4.1.1             lazyeval_0.2.2            spatstat.geom_2.2-2       yaml_2.2.1               
#[56] reshape2_1.4.4            abind_1.4-5               httpuv_1.6.3              tools_4.1.1               ellipsis_0.3.2            spatstat.core_2.3-0       RColorBrewer_1.1-2        ggridges_0.5.3            Rcpp_1.0.7                plyr_1.8.6                sparseMatrixStats_1.4.2  
#[67] zlibbioc_1.38.0           purrr_0.3.4               RCurl_1.98-1.4            rpart_4.1-15              deldir_0.2-10             pbapply_1.5-0             viridis_0.6.1             cowplot_1.1.1             zoo_1.8-9                 SeuratObject_4.0.2        ggrepel_0.9.1            
#[78] cluster_2.1.2             magrittr_2.0.1            data.table_1.14.0         scattermore_0.7           lmtest_0.9-38             RANN_2.6.1                fitdistrplus_1.1-5        patchwork_1.1.1           mime_0.11                 evaluate_0.14             xtable_1.8-4             
#[89] gridExtra_2.3             compiler_4.1.1            tibble_3.1.4              KernSmooth_2.23-20        crayon_1.4.1              R.oo_1.24.0               htmltools_0.5.2           mgcv_1.8-36               later_1.3.0               tidyr_1.1.3               DBI_1.1.1                
#[100] tweenr_1.0.2              dbplyr_2.1.1              MASS_7.3-54               Matrix_1.3-4              permute_0.9-5             cli_3.0.1                 R.methodsS3_1.8.1         igraph_1.2.6              pkgconfig_2.0.3           plotly_4.9.4.1            scuttle_1.2.1            
#[111] spatstat.sparse_2.0-0     vipor_0.4.5               dqrng_0.3.0               XVector_0.32.0            stringr_1.4.0             digest_0.6.27             sctransform_0.3.2         RcppAnnoy_0.0.19          vegan_2.5-7               spatstat.data_2.1-0       rmarkdown_2.11           
#[122] leiden_0.3.9              uwot_0.1.10               edgeR_3.34.1              DelayedMatrixStats_1.14.3 shiny_1.6.0               gtools_3.9.2              lifecycle_1.0.0           nlme_3.1-152              jsonlite_1.7.2            Rhdf5lib_1.14.2           BiocNeighbors_1.10.0     
#[133] viridisLite_0.4.0         limma_3.48.3              fansi_0.5.0               pillar_1.6.2              lattice_0.20-44           fastmap_1.1.0             httr_1.4.2                survival_3.2-13           glue_1.4.2                png_0.1-7                 ggforce_0.3.3            
#[144] stringi_1.7.4             HDF5Array_1.20.0          BiocSingular_1.8.1        dplyr_1.0.7               irlba_2.3.3               future.apply_1.8.1      