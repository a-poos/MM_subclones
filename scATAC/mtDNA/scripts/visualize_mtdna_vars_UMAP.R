#############################################################################################################################
##                                                                                                                      
##  VISUALIZE MTDNA VARIANTS ON THE ARCHR UMAPS
##                                                                                                                      
##  Date: 02 APRIL 2021, modified 06 July 2022                                                                                                              
##  
##  Authors: Moritz Przybilla, Alexandra Poos
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("devtools", "BiocManager", "reshape2", "ArchR", "tidyverse", "Matrix", "harmony", "pheatmap", "cowplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors) but check whether everything is in place
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# FUNCTIONS
#####################################################################################
"%ni%" <- Negate("%in%")

source("/home/bq_apoos/sd19k002/scATAC-seq/Github/scATAC/mtDNA/scripts/variant_calling.R")

plot_mito_umap <- function(mito_df, variant){
  
  # Ugly name that comes from making it a column names
  ugly_var <- variant
  mito_df$variant_value <- mito_df[,ugly_var]*100
  
  ggplot(mito_df %>% arrange(variant_value), aes(x = UMAP1, y = UMAP2, color = variant_value)) +
    geom_point() + theme_classic() +
    xlab("") + 
    ylab("") +
    scale_color_gradientn(colors = colorRampPalette(c("lightgrey","firebrick4"))(50)) +
    labs(color = variant) +
    theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
  
}

############################################################################
##                          READ IN THE DATA
############################################################################

# define in and output directory
input.dir <- "/home/bq_apoos/sd19k002/scATAC-seq/Github/scATAC/mtDNA/mgatk"
output.dir <- "/home/bq_apoos/sd19k002/scATAC-seq/Github/scATAC/mtDNA/plots"
dir.create(paste0(output.dir, "/UMAP_mtdna"))

# list the archr projects
archr.projects <- list.files("/home/bq_apoos/sd19k002/scATAC-seq/Github/scATAC/ArchR", full.names = T)

# mgatk_se is the Summarized Experiment .rds file
# This is automatically produced from running mgatk.
se.files <- list.files(input.dir, pattern = ".rds", recursive = T, full.names = T, all.files = T)
#se.files <- se.files[grep("T1/|T2/", se.files)]

# get sample ids
sample.ids <- str_split_fixed(basename(se.files), ".rds", 2)[,1]

# get the patient ids
patient.ids <- unique(str_split_fixed(sample.ids, "-", 2)[,1])

# read in the filtered mutation dataframe
mgatk.variants <- read.table("/home/bq_apoos/sd19k002/scATAC-seq/Github/scATAC/mtdna/mgatk/mtDNA_variants_RRMM15.txt", header = T)

# filter low QC variants
mgatk.variants <- mgatk.variants[round(mgatk.variants$mean, 2) >= 0.01,]

i <- 1
for (i in 1:length(patient.ids)){
  
  # read in one of the mgatk SE objects
  patient.tmp <- patient.ids[i]
  cat(paste0("\nProcessing patient...", patient.tmp, "\n"))
  samples.tmp <- sample.ids[grep(patient.tmp, sample.ids)]
  
  if(length(samples.tmp) > 1){
    
    pt.se.files <- se.files[grep(patient.tmp, se.files)]
    mgatk_se <- lapply(pt.se.files, readRDS)
    
    #colnames(mgatk_se[[1]])=paste(patient.ids[i], "-T1#", colnames(mgatk_se[[1]]), sep="")
    #colnames(mgatk_se[[2]])=paste(patient.ids[i], "-T2#", colnames(mgatk_se[[2]]), sep="")
    
    #for RRMM15 as example
    colnames(mgatk_se[[1]])=paste("TS6ZX9-T1#", colnames(mgatk_se[[1]]), sep="")
    colnames(mgatk_se[[2]])=paste("TS6ZX9-T2#", colnames(mgatk_se[[2]]), sep="")
    
    mgatk_se <- cbind(mgatk_se[[1]], mgatk_se[[2]])
    
  } else {
    
    pt.se.files <- se.files[grep(patient.tmp, se.files)]
    mgatk_se <- readRDS(pt.se.files)
    colnames(mgatk_se)=paste(patient.ids[i], "-T2#", colnames(mgatk_se), sep="")
  }
  
  # call mutations in the mtDNA
  mut_se <- call_mutations_mgatk(mgatk_se)
  
  # grep the patient variants and filter matrix
  pt.vars <- unique(as.character(mgatk.variants[mgatk.variants$DonorID == patient.tmp, "variant"]))
  
  # LOAD THE PROJECT OF INTEREST
  projRRMM <- loadArchRProject(archr.projects[grep(patient.tmp, archr.projects)])
  
  # get the umap coordinates
  umap.df <- projRRMM@embeddings$UMAP$df
  colnames(umap.df) <- c("UMAP1", "UMAP2")
  umap.df$Cell_barcodes=rownames(umap.df)
  
  # get the vaf and coverage per cell per mutation
  vaf.df <- as.data.frame(as.matrix(t(assays(mut_se)[["allele_frequency"]][as.character(pt.vars),])))
  cov.df <- as.data.frame(as.matrix(t(assays(mut_se)[["coverage"]][as.character(pt.vars),])))
  colnames(cov.df) <- paste0(colnames(cov.df), "_cov")
  
  # combine both
  combined.df <- cbind(vaf.df, cov.df)
  combined.df$Cell_barcodes=rownames(combined.df)
  
  # combine to mito dataframe
  mito_df <- merge(umap.df, combined.df, by = "Cell_barcodes")

  ############################################################################
  ##                          VISUALIZE UMAP
  ############################################################################
  
  # plot UMAP by clusters and samples
  # color by Sample
  p1 <- plotEmbedding(ArchRProj = projRRMM, colorBy = "cellColData", name = "sample", embedding = "UMAP", size = 2)
  
  # color by Clusters
  p2 <- plotEmbedding(ArchRProj = projRRMM, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 2)
  
  # align the individual plots
  pdf(paste0(output.dir, "/UMAP_mtdna/", patient.tmp, "_UMAP_sample_clusters.pdf"), width = 10, height = 10)
  print(ggAlignPlots(p1, p2, type = "h"))
  dev.off()
  
  bla=colnames(mito_df)[grep("_cov", colnames(mito_df))]
  bla=gsub("_cov", "", bla)
  mito_df2=mito_df
  mito_df2[,bla][mito_df2[,bla]>0.1] <- 0.1
  
  j <- 1
  for (j in 1:length(pt.vars)){
    
    variant <- pt.vars[j]
    plot <- plot_mito_umap(mito_df2, variant)
    #plot <- plot_mito_umap(mito_df, variant)
  
    pdf(paste0(output.dir, "/UMAP_mtdna/", patient.tmp, "_UMAP_", variant, ".pdf"), width = 6, height = 5)
    print(plot)
    dev.off()
    
    projRRMM@cellColData$mvariant="NA"
	for (ii in 1:length(rownames(projRRMM@cellColData))){
    	if(rownames(projRRMM@cellColData)[ii] %in% mito_df$Cell_barcode){
        	x=which(mito_df$Cell_barcode==rownames(projRRMM@cellColData)[ii])
        	projRRMM@cellColData$mvariant[ii]=as.numeric(mito_df[x,which(colnames(mito_df)==pt.vars[j])])
    	}
	}
	
    colnames(projRRMM@cellColData)[which(colnames(projRRMM@cellColData)=="mvariant")]=paste("Mito", (gsub(">", "", pt.vars[j])), sep="_")
  }
}


#### Heatmap
library(BuenColors)
library(circlize)
library(ComplexHeatmap)
library(SummarizedExperiment)


# exemplarily for RRMM15
rrmm15=cbind(projRRMM$cellNames, projRRMM$Mito_10114TC, projRRMM$Mito_12962GA, projRRMM$Mito_13681AG, projRRMM$Mito_13958GA, projRRMM$Mito_14285TC, projRRMM$Mito_14769AG, projRRMM$Mito_15356GA,  projRRMM$Mito_15928GA, projRRMM$Mito_5115TC, projRRMM$Mito_821TC, projRRMM$Mito_9868GA, projRRMM$clones_int)

 colnames(rrmm15)=c("cellNames", "Mito_10114TC", "Mito_12962GA", "Mito_13681AG", "Mito_13958GA", "Mito_14285TC", "Mito_14769AG", "Mito_15356GA",  "Mito_15928GA", "Mito_5115TC", "Mito_821TC", "Mito_9868GA", "Subclone")
 rownames(rrmm15)=rrmm15[,1]
 rrmm15=rrmm15[,-1]
 
 rrmm15=rrmm15[order(rrmm15[,12], decreasing=F),]
 
 test=rrmm15
 test=test[order(test[,12], decreasing=F),]
 test=test[,-12]
 test=as.matrix(test)
 test=t(test)
 class(test)<-"numeric"
 
 test2=test
 test2[test2>0.1] <- 0.1
 
 color_vec <- c("#D51F26", "#F47D2B", "#272E6A", "#89288F", "#8A9FD1", "#C06CAB", "#D8A767", "#FEE500","#90D5E4", "#89C75F", "#F37B7D"); names(color_vec) <- c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5", "Clone6", "Clone7", "Clone8", "Clone9", "Clone10", "Clone11")
 
ha_col <- HeatmapAnnotation(df = data.frame(Subclone =rrmm15[,12]),
                           col = list(Subclone = color_vec)) 
 
hm <- Heatmap(data.matrix(test), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              top_annotation=ha_col,
              name = "AF",
              cluster_rows = FALSE, 
              show_column_names = FALSE)


hm2 <- Heatmap(data.matrix(test2), 
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              cluster_columns = FALSE,
              top_annotation=ha_col,
              name = "AF",
              cluster_rows = FALSE, 
              show_column_names = FALSE)


pdf(paste0(output.dir, "/UMAP_mtdna/", "Heatmap_Mitos_per_Cell_RRMM15.pdf"))
hm
dev.off()

pdf(paste0(output.dir, "/UMAP_mtdna/", "Heatmap_Mitos_per_Cell_RRMM15_01.pdf"))
hm2
dev.off()

