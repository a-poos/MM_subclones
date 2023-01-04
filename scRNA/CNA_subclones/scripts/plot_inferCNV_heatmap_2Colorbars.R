#############################################################################################################################
##                                                                                                                      
##  VISUALIZE THE INFERCNV COPY NUMBER HEATMAP
##                                                                                                                      
##  Date: 12 NOVEMBER 2020, modified 16 December 2021                                                                                                          
##  
##  Author: Moritz Przybilla and Alexandra Poos
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", 
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "BiocManager", 
                      "biomaRt", "httr", "data.table", "Seurat", "stringr", "fields", "dplyr", "GenomicRanges", "phylogram", "dendextend")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
httr::set_config(httr::config(ssl_verifypeer=0L))

############################################################################
##                              FUNCTIONS
############################################################################

# negate %in%
"%ni%" <- Negate("%in%")

# GENE EXPRESSION WITH TIME POINT MATRIX PLOT
plot_geneExp_time_heatmap <- function(norm.gExp.matrix, gene.ordering.file, cluster_dataframe, sample_name, output_directory = "/home/bq_apoos/sd19k002/scRNA-seq/PairedSamples/Analysis_CR5/InferCNV/Biopsy_PCs_C14_single/Heatmaps/"){
  
  # get the clusters for the respective sample stored as a dataframe with a clone_id column
  clusters <- as.numeric(str_split_fixed(unique(cluster_dataframe$clone_id), "Clone", 2)[,2])
  valid.clones <- factor(clusters, levels = c(1:as.numeric(max(unique(clusters)))))
  valid.clones <- valid.clones[order(valid.clones)]
  
  # center the infercnv output to 0
  # gexp.norm <- norm.gExp.matrix -1
  gexp.norm <- norm.gExp.matrix
  
  # only take genes which are present in the matrix
  rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol
  genes <- makeGRangesFromDataFrame(gene.ordering.file)
  genes <- genes[names(genes) %in% rownames(gexp.norm)]
  
  # generate new object genes.interest with genes and respective coordinates
  genes.interest <- as.data.frame(genes)
  rownames(genes.interest) <- names(genes)
  mat <- gexp.norm
  
  # check for each chromosome if there are genes present
  gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
    na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  chrs <- paste0("chr", 1:22)
  gExp.array <- gExp.array[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100
  
  # generate the layout accordingly
  plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)
  
  # here we define the groups we want to cluster according to
  groups <- valid.clones
  chr.list <- list()
  ordered_names <- c()
  j <- 1
  
  # perform hierarchical clustering within each of these groups
  for (j in 1:length(groups)){
    
    # which group?
    print(groups[j])
    
    # get the barcodes from the respective group
    barcodes <- cluster_dataframe[grep(paste("Clone",groups[j],"$", sep=""), cluster_dataframe$clone_id), "Cell_barcodes"]
    barcodes <- barcodes$Cell_barcodes
    
    if (length(barcodes) > 1){
      # iterate over all chrs and subset the matrix
      
      for (i in 1:length(gExp.array)){
        
        # matrix 
        matrix <- gExp.array[[i]]
        matrix <- matrix[,colnames(matrix) %in% barcodes]
        
        chr.list[[paste0("chr",i)]] <- matrix
        
      }
      
      # calculate the col averages here for the subset
      avgd <- do.call(rbind, lapply(names(chr.list),function(nam) {
        d <- chr.list[[nam]]
        d <- colMeans(d)
        d
      }))
      
      # Order cells with hierarchical clustering
      dist.centered.matrix <- dist(as.matrix(t(avgd)), method = "euclidean")
      hc <- hclust(dist.centered.matrix, method = "ward.D2")
      
      # make a vector with the right order of barcodes per group
      ordered_names <- c(ordered_names, hc$labels[hc$order]) 
    } else {
      
      # add single barcode to the list
      ordered_names <- c(ordered_names, barcodes)
      cat(paste0("Clone", groups[j], " does only have 1 cell!\n"))
      next
    }
    
  }
  
  # set parameters for the image
  adapted.gExp.list <- gExp.array
  pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
  zlim <- c(0.9, 1.1)
  
  # limit the gExp to maximums according to zlim
  limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
    d <- adapted.gExp.list[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  #####################################################################################
  #                   GENERATE THE COLOUR BAR FOR THE HEATMAP
  #####################################################################################
  
  # first, remove the clone id and create and index
  cluster_dataframe$clone_number <- as.numeric(str_split_fixed(cluster_dataframe$clone_id, "Clone", 2)[,2])
  
  # the index is important for the order so we merge it with the seurat metadata
  #heatmap.metadata <- cluster_dataframe[,c("Cell_barcodes", "clone_number", "timepoint")]
  heatmap.metadata <- cluster_dataframe[,c("Cell_barcodes", "clone_number", "group")]
  rownames(heatmap.metadata) <- heatmap.metadata$Cell_barcodes
  heatmap.metadata <- heatmap.metadata[order(heatmap.metadata$group),]
  heatmap.metadata <- heatmap.metadata[order(heatmap.metadata$clone_number),]
  
  # add a colour column for the time points
  heatmap.metadata$colour <- "colour"
  heatmap.metadata[grep("T1", heatmap.metadata$group),"colour"] <- "grey80" #blue
  heatmap.metadata[grep("T2", heatmap.metadata$group),"colour"] <- "grey60" #red  
  
  # add a col name and a index for the cells to monitor each cells position in the matrix
  rownames(heatmap.metadata) <- heatmap.metadata$Cell_barcodes
  heatmap.metadata$final_index <- c(1:nrow(heatmap.metadata))
  heatmap.metadata[grep("T1", heatmap.metadata$group), "final_index"] <- 1
  heatmap.metadata[grep("T2", heatmap.metadata$group), "final_index"] <- 2
  
  # make a df with the index to color it accordingly
  heatmap.timepoint.col.bar <- data.frame(as.numeric(heatmap.metadata$final_index),as.numeric(heatmap.metadata$final_index))
  heatmap.clone.col.bar <- data.frame(as.numeric(heatmap.metadata$clone_number),as.numeric(heatmap.metadata$clone_number))
  colnames(heatmap.clone.col.bar)<-c("clone1","clone2")
  
  # set the colors
  colors <- c("#D51F26", "#F47D2B", "#272E6A", "#89288F", "#8A9FD1", "#C06CAB", "#D8A767", "#FEE500","#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8")
  
   # get the coloring according to RGs
  cols.timepoint <- c(as.character(unique(heatmap.metadata$colour)))
  cols.clones <- colors[1:length(valid.clones)]
  
  #####################################################################################
  #                 GENERATE THE HEATMAP ITSELF AND WRITE TO FILE
  #####################################################################################

  pdf(paste0(output_directory, "/", sample_name, "/", sample_name, "_inferCNV_timepoint_heatmap.pdf") , width = 26, height = 12)
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  image(t(heatmap.clone.col.bar),xlab="",ylab="", axes=FALSE, col=cols.clones)
  box()
  image(t(heatmap.timepoint.col.bar),xlab="",ylab="", axes=FALSE, col=cols.timepoint)
  ## plot chromosomes
  box()
  for (i in 1:length(limit.gExp.list)){
    message(chrs[i])
    d <- limit.gExp.list[[i]]
    d <- d[, heatmap.metadata$Cell_barcodes] 
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=chrs[i], useRaster = T)
    box()
  }
  dev.off()
  
  png(paste0(output_directory, "/", sample_name, "/", sample_name, "_inferCNV_heatmap_legend.png") , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
  par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
  image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
  dev.off()
  
}

############################################################################
##                  READ IN ALL THE FILES OF INTEREST
############################################################################
# set sample id
sample.tmp <- "RRMM15"

# specify the input directory
input.dir <- "/Github/scRNA/CNA_subclones/input"
# specify the output directory
out.dir <- "/Github/scRNA/CNA_subclones/plots"

# create a sample directory
dir.create(paste0(out.dir, "/", sample.tmp))

## GET ALL THE SAMPLE INFORMATION REQUIRED FOR 10X RNA DATA FROM INFERCNV
# list all scRNA matrices
scRNA.files <- list.files(input.dir, pattern = "infercnv.observations.txt", full.names = T, recursive = T)
scRNA.file <- scRNA.files[grep(sample.tmp, scRNA.files)]

# get the scRNA clusters files
scRNA.cluster.file <- paste0(input.dir, "/","RRMM15_11subclones.txt")

############################################################################
##                  VISUALIZE THE INFERCNV HEATMAP
############################################################################

# read in the scRNA matrix
norm.gExp.matrix = as.matrix(read.table(scRNA.file, header = T, sep = " "))
norm.gExp.matrix=norm.gExp.matrix[,-grep("_N", colnames(norm.gExp.matrix))] #remove all cells from the CD138 negative fraction

# and gene ordering files
gene.order.file <- paste0(input.dir, "/","gene_ordering_file.txt")
gene.order.file <- read.table(gene.order.file, header = F, col.names = c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id"))
gene.ordering.file = gene.order.file

cluster_dataframe = fread(scRNA.cluster.file, header = T, sep = "\t", col.names = c("Cell_barcodes", "clone_id", "group"))
sample_name = sample.tmp 

# plot CNA heatmap
plot_geneExp_time_heatmap(norm.gExp.matrix, gene.ordering.file, cluster_dataframe, sample_name,  output_directory = out.dir)


############################################################################
##                  EXAMPLE of the WGS-guided clustering approach
############################################################################
# cluster cells of RRMM15 T2 based on chromosome3
norm.gExp.matrix=norm.gExp.matrix[,grep("_T2_", colnames(norm.gExp.matrix))]
selected=gene.order.file[gene.order.file$chr %in% c("chr3"),]
#### if necessary specify specific region by selected=selected[selected$start>75000000,] or selected=selected[selected$start<75000000,]

# cluster cells only based on chr3
norm.gExp.matrix3=norm.gExp.matrix
norm.gExp.matrix3=norm.gExp.matrix3[rownames(norm.gExp.matrix3) %in% selected$hgnc_symbol,]

p=hclust(dist(t(norm.gExp.matrix3)), method="ward.D2")
plot(p,labels = FALSE,hang = -1,xlab = "scRNA",ylab = "distance")

NCLONES <- 3

tp <- color_branches(p, k = NCLONES)
plot(tp,leaflab = 'none')

g <- cutree(p,k = NCLONES)
names(g)=gsub("_2", "", names(g))
groups <- data.frame(BARCODE=names(g),#BARCODE=barcodes,
                     CLONE=as.numeric(g),
                     stringsAsFactors = F)

colnames(groups)=c("Cell_barcodes", "clone_id")
groups$clone_id <- paste0("Clone",groups$clone_id)

groups$group="T1"
groups$group[grep("T2", groups$Cell_barcodes)]="T2"


write.table(groups, paste0(input.dir, "/","RRMM15_T2_3subclones_chr3.txt"), sep="\t", col.names=T, row.names=F, quote=F)
groups=fread(paste0(input.dir, "/","RRMM15_T2_3subclones_chr3.txt"), header = T, sep = "\t", col.names = c("Cell_barcodes", "clone_id", "group"))


plot_geneExp_time_heatmap(norm.gExp.matrix, gene.ordering.file, groups, sample_name,  output_directory = out.dir)
   