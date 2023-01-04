#############################################################################################################################
##                                                                                                                      
##  VISUALIZE THE ATAC COPY NUMBER HEATMAP
##                                                                                                                      
##  Date: 17 November 2021                                                                                                          
##  
##  Author: Moritz Przybilla and Alexandra Poos
##
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", "SummarizedExperiment", "BuenColors", "circlize",
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "BiocManager", "ggbeeswarm", "ComplexHeatmap",
                      "biomaRt", "httr", "data.table", "Seurat", "stringr", "fields", "dplyr", "GenomicRanges", "phylogram", "dendextend")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
httr::set_config(httr::config(ssl_verifypeer=0L))

# negate %in%
"%ni%" <- Negate("%in%")


###############################################################################################
##       calculate CNA profiles based on Lareau et al. (Nature Biotechnology, 2021)
###############################################################################################

# Get the typical profile of a fragment distribution
normal_npcs <- readRDS("Normal_PCs_cnv.rds"); normal_npcs[is.na(normal_npcs)] <- 0

cpm_norm <- (t(t(normal_npcs)/colSums(normal_npcs)) * 100000)
row_means <- rowMeans(cpm_norm)
row_std <- sqrt(rowVars(cpm_norm))

mat <- readRDS("ATAC_TS6ZX9_cnv.rds"); mat[is.na(mat)] <- 0
makeZscoreMat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- (mat_cpm_norm - row_means)/row_std
}

makeLog2Mat <- function(mat){
  mat_cpm_norm <- (t(t(mat)/colSums(mat)) * 100000)
  zscore_mat <- log2(mat_cpm_norm/(row_medians + 1))
  zscore_mat
}

zscore <- makeZscoreMat(mat)
score_bs <- makeZscoreMat(normal_npcs)
region_meta <- str_split_fixed(rownames(zscore), "_", 4)

keep <- TRUE
ordering <- factor(region_meta[keep,1], levels = unique(region_meta[keep,1]))

zscore[zscore > 3] <- 3
zscore[zscore < -3] <- -3


###############################################################################################
##       			visualize CNA heatmap, RRMM15 as example
###############################################################################################

  groups2=fread("RRMM15_11subclones.txt", header = T, sep = "\t")
  groups2$Cell_barcodes <- str_split_fixed(groups2$Cell_barcode, "#", 2)[,2]

  mito_df=groups2
  mito_df$mito_clone13="T1"
  mito_df$mito_clone13[grep("T2", mito_df$Cell_barcode)]="T2"

  cluster_dataframe=mito_df
  mito="mito_clone13"
  sample_name="RRMM15"
  output_directory="/home/bq_apoos/sd19k002/scATAC-seq/Github/scATAC/CNA_subclones/plots"


  cluster_dataframe$clone_number <- as.numeric(str_split_fixed(cluster_dataframe$clone_id, "Clone", 2)[,2])
  
  # get the clusters for the respective sample stored as a dataframe with a clone_id column
  clusters <- as.numeric(str_split_fixed(unique(cluster_dataframe$clone_id), "Clone", 2)[,2])
  valid.clones <- factor(clusters, levels = c(1:as.numeric(max(unique(clusters)))))
  valid.clones <- valid.clones[order(valid.clones)]
  valid.clones <- valid.clones[!is.na(valid.clones)]
  
  # convert cluster dataframe
  cluster_dataframe <- as.data.frame(cluster_dataframe)
  
  # center the infercnv output to 0
  atac.mtx <- zscore
  
  # generate new object genes.interest with genes and respective coordinates
  mat <- as.data.frame(atac.mtx)
  
  # check for each chromosome if there are genes present  
  test=list(zscore[grep("chr1_", rownames(zscore)),], zscore[grep("chr2_", rownames(zscore)),], zscore[grep("chr3_", rownames(zscore)),], zscore[grep("chr4_", rownames(zscore)),], zscore[grep("chr5_", rownames(zscore)),], zscore[grep("chr6_", rownames(zscore)),], zscore[grep("chr7_", rownames(zscore)),], zscore[grep("chr8_", rownames(zscore)),], zscore[grep("chr9_", rownames(zscore)),], zscore[grep("chr10_", rownames(zscore)),], zscore[grep("chr11_", rownames(zscore)),], zscore[grep("chr12_", rownames(zscore)),], zscore[grep("chr13_", rownames(zscore)),], zscore[grep("chr14_", rownames(zscore)),], zscore[grep("chr15_", rownames(zscore)),], zscore[grep("chr16_", rownames(zscore)),], zscore[grep("chr17_", rownames(zscore)),], zscore[grep("chr18_", rownames(zscore)),], zscore[grep("chr19_", rownames(zscore)),], zscore[grep("chr20_", rownames(zscore)),], zscore[grep("chr21_", rownames(zscore)),], zscore[grep("chr22_", rownames(zscore)),])
  chrs <- paste0("chr", 1:22)
  names(test)=chrs
  
  gatac.array=test
  
  # we want to scale the output to the respective amount of genes and their widths
  widths <- sapply(gatac.array, nrow); widths <- widths/max(widths)*100
  
  # generate the layout accordingly
  plot.layout <- layout(matrix(seq(1,length(gatac.array)),1,length(gatac.array),byrow=TRUE), widths=widths)
  
  # here we define the groups we want to cluster according to
  groups <- c("T1", "T2")
  clones <- valid.clones
  chr.list <- list()
  ordered_names <- c()
  group.df <- cluster_dataframe[,c(1,3, grep("clone_id", colnames(cluster_dataframe)),grep("mito_clone13", colnames(cluster_dataframe)))]
  j <- 1
  
  # perform hierarchical clustering within each of these groups
  for (j in 1:length(clones)){
    
    # which group
    print(clones[j])
    
    # get the barcodes from the respective group
    #group.mito.df <- group.df[grep(clones[j], group.df[,"clone_id"]),]
    group.mito.df <- group.df[group.df$clone_id %in% paste0("Clone", clones[j]),]
    print(nrow(group.mito.df))
    clone_ordered_names <- c()
    k <- 1
    
    for (k in 1:length(groups)){
      
      #print(clones[k])
      print(groups[k])
      
      # get the barcodes from the respective subgroup
      barcodes <- group.mito.df[grep(groups[k], group.mito.df$Cell_barcode),"Cell_barcode"]
      print(length(barcodes))
      
      if (length(barcodes) > 1){
        # iterate over all chrs and subset the matrix
        for (i in 1:length(gatac.array)){
          
          # matrix 
          matrix <- gatac.array[[i]]
          matrix=matrix[,colnames(matrix) %in% barcodes]
          
          chr.list[[paste0("chr",i)]] <- matrix
        }
        
        # calculate the col averages here for the subset
        avgd <- do.call(rbind, lapply(names(chr.list),function(nam) {
          d <- chr.list[[nam]]
          d <- colMeans(d)
          d
        }))
        
        # Order cells with hierarchical clustering
        dist.centered.matrix <- dist(t(avgd), method = "euclidean")
        hc <- hclust(dist.centered.matrix, method = "ward.D2")
        
        # make a vector with the right order of barcodes per group
        clone_ordered_names <- c(clone_ordered_names, hc$labels[hc$order])
      }
    }
    
    # make a vector with the right order of barcodes per group
    ordered_names <- c(ordered_names, clone_ordered_names)
  }
  
  # set parameters for the image
  adapted.gatac.list <- gatac.array
  pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
  # pcol_corr <- colorRampPalette(c("white", "darkgreen"))(256)
  zlim <- c(-3, 3) # upper and lower limit
  # zlim_corr <- c(0, max(cluster_dataframe$pearson.correlation))
  
  # limit the gExp to maximums according to zlim
  limit.gatac.list <- lapply(names(adapted.gatac.list),function(nam) {
    d <- adapted.gatac.list[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  #####################################################################################
  #                   GENERATE THE COLOUR BAR FOR THE HEATMAP
  #####################################################################################
  
  # the index is important for the order so we merge it with the seurat metadata
  clone.mito.annotation <- cluster_dataframe[cluster_dataframe$Cell_barcode %in% ordered_names,]
  clone.mito.annotation <- clone.mito.annotation[,c("Cell_barcode", "clone_number", mito)]
  rownames(clone.mito.annotation) <- clone.mito.annotation$Cell_barcode
  clone.mito.annotation <- clone.mito.annotation[ordered_names,]
  
  # add a col name and a index for the cells to monitor each cells position in the matrix
  clone.mito.annotation$final_index <- c(1:nrow(clone.mito.annotation))
  clone.mito.annotation[grep("T1", clone.mito.annotation[,mito]), "final_index"] <- 1
  clone.mito.annotation[grep("T2", clone.mito.annotation[,mito]), "final_index"] <- 2
  
  # assign sample id
  clone.mito.annotation$colour <- "colour"
  clone.mito.annotation[grep("T1", clone.mito.annotation[,mito]),"colour"] <- "grey80"
  clone.mito.annotation[grep("T2", clone.mito.annotation[,mito]),"colour"] <- "grey60"
  
  # make a df with the index to color it accordingly
  mito.clone.col.bar <- data.frame(as.numeric(clone.mito.annotation$final_index),as.numeric(clone.mito.annotation$final_index))
  cnv.clone.col.bar <- data.frame(as.numeric(clone.mito.annotation$clone_number),as.numeric(clone.mito.annotation$clone_number))
  
  # set the colors
 colors <- c("#D51F26", "#F47D2B", "#272E6A", "#89288F", "#8A9FD1", "#C06CAB", "#D8A767", "#FEE500","#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8")
  
  mito.cols <- unique(clone.mito.annotation$colour)
  cnv.clone.cols <- colors[1:length(unique(valid.clones))]
  
  #####################################################################################
  #                 GENERATE THE HEATMAP ITSELF AND WRITE TO FILE
  #####################################################################################
  
  pdf(paste0(output_directory, "/", sample_name, "/", sample_name, "_scATAC_CNAsubclones_heatmap.pdf") , width = 26, height = 12)
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(t(cnv.clone.col.bar),xlab="",ylab="", axes=FALSE, col=cnv.clone.cols)
  box()
  image(t(mito.clone.col.bar),xlab="",ylab="", axes=FALSE, col=mito.cols)
  ## plot chromosomes
  box()
  for (i in 1:length(limit.gatac.list)){
    message(chrs[i])
    d <- limit.gatac.list[[i]]
    d <- as.matrix(d[, clone.mito.annotation$Cell_barcode])
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=chrs[i], useRaster = T)
    box()
  }
  # box()
  dev.off()
  
  png(paste0(output_directory, "/", sample_name, "/", sample_name, "_scATAC_CNAsubclones_heatmap_legend.png") , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
  image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
  dev.off()
 
 
 
 

###############################################################################################
##       			WGS guided clustering approach, RRMM15 T2 as example
###############################################################################################  

c3=region_meta[which(region_meta[,1]=="chr3"),]

selected=c3
selected=as.data.frame(selected)
selected$V5=paste(selected$V1, selected$V2, selected$V3, selected$V4, sep="_")
zscore2=zscore[rownames(zscore) %in% selected$V5,]
zscore2=zscore2[,colnames(zscore2) %in% cl24$Cell_barcode]

p=hclust(dist(t(zscore2)), method="ward.D2")
plot(p,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

NCLONES <- 3

tp <- color_branches(p, k = NCLONES)
plot(tp,leaflab = 'none')

g <- cutree(p,k = NCLONES)
names(g)=gsub("_2", "", names(g))
groups <- data.frame(BARCODE=names(g),#BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)

colnames(groups)=c("Cell_barcode", "clone_id")
groups$clone_id <- paste0("Clone",groups$clone_id)

write.table(groups, "RRMM15-T2_3subclones_chr3.txt", sep="\t", col.names=T, row.names=F, quote=F)

groups2=fread("RRMM15-T2_3subclones_chr3.txt", header = T, sep = "\t", col.names = c("Cell_barcode", "clone_id"))
groups2$Cell_barcodes <- str_split_fixed(groups2$Cell_barcode, "#", 2)[,2]

mito_df=groups2
mito_df$mito_clone13="T1"
mito_df$mito_clone13[grep("T2", mito_df$Cell_barcode)]="T2"

cluster_dataframe=mito_df


#### visualize heatmap again (see section above)
