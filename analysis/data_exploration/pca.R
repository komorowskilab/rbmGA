#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

set_wd_to_script_dir <- function() {
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && !is.null(rstudioapi::getActiveDocumentContext()$path) &&
      rstudioapi::getActiveDocumentContext()$path != "") {
    path <- rstudioapi::getActiveDocumentContext()$path
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_flag <- "--file="
    match <- grep(file_flag, args)
    if (length(match) > 0) {
      path <- sub(file_flag, "", args[match])
    } else {
      stop("Could not determine script path. Run via Rscript or from an RStudio-opened file.")
    }
  }
  setwd(dirname(normalizePath(path)))
  message("Working directory set to: ", getwd())
}

set_wd_to_script_dir()
getwd()
set.seed(0)

#load librarires
if (!require(data.table)) install.packages('data.table')
library(data.table)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

timepoint='D7'

#default is set to 500 for this particular project
intersecting_colnames<-function(val=500){
  #MCFS selection
  res<-readRDS('../../data/Mcfs_results/mcfs_D7_gp1.rds')
  res_1<-readRDS('../../data/Mcfs_results/mcfs_D133_gp1.rds')
  sig_genes<- res$RI$attribute[1:val]
  sig_genes_1<- res_1$RI$attribute[1:val]
  union(sig_genes,sig_genes_1)
  intersect(sig_genes,sig_genes_1)
  col<-intersect(sig_genes,sig_genes_1)
  return(col)
}

read_data<-function(x,col){
  print(x)
  #reading the Day 7 data
  data<-fread(paste0('../../data/Rosetta_Decision_Tables/',x,'.csv'))

  #moving animal ids to rownames inorder to keep information intact
  data <- data %>% column_to_rownames(.,'AnimalID')

  #data preparation by subsetting the data with 70 genes and outcome
  df<-subset(data, select=c(col,'Group','protectionStatus'))
  df<-df[df$Group =='S' | df$Group =='O' | df$Group =='ABL1_a',]
  ros_data<- subset(df, select=c(col,'protectionStatus')) #

  return(ros_data)
}


#Function to execute PCA; pcadf is the df underinvestigation, mcfs_res isthe mcfs res for the binary classificatio
# subset_col to be used to invesigate
pca_fun<-function(pcadf,mcfs_res,subset_col=FALSE){
  library(ggfortify)
  library("FactoMineR")
  library("factoextra")
  if (subset_col==TRUE){
    cols<- mcfs_res$RI[1:mcfs_res$cutoff_value,'attribute']
    colnames(pcadf)<- gsub('-','_',colnames(pcadf))
    pcadf<- subset(pcadf,select=c(cols,'attr_oncotree_disease_code'))

  }
  df.pca <- prcomp(pcadf[,2:(ncol(pcadf)-1)],center = TRUE,scale. = TRUE)
  #class<-as.factor(pcadf$)

  #Diferent PCA plots
  # pca.plot <- autoplot(df.pca,data = pcadf,colour = 'protectionStatus')+ ggtitle(ncol(pcadf))
  # print(pca.plot)
  # print(fviz_pca_ind(df.pca, col.ind = "contrib", xlab='PC1',ylab='PC2',
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #              repel = TRUE # Avoid text overlapping (slow if many points)
  # ))
  # print(fviz_pca_var(df.pca, col.var="contrib",label='var',xlab='PC1',ylab='PC2',
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #              repel = TRUE, # Avoid text overlapping
  #              axes = c(1, 2) ))# choose PCs to plot
  #
  # print(fviz_pca_biplot(df.pca, repel = TRUE,
  #                       col.var = "#2E9FDF", # Variables color
  #                       col.ind = "#696969"  # Individuals color
  # ))

  #for manuscript
  print(fviz_eig(df.pca, addlabels = TRUE, ylim = c(0, 50),ncp=15))


  print(fviz_pca_ind(df.pca,
                        # Individuals
                        geom.ind = "point",
                        fill.ind = pcadf[,ncol(pcadf)],
                         palette = c("#999999", "#E69F00"),title='',
                       # palette = "jco",
                        addEllipses = TRUE, col.var = "contrib",ellipse.type = "convex",
                        legend.title = list(fill = "Protection-Status" )))


  manuscript_plot_pca<- fviz_pca_ind(df.pca, geom.ind = "point",
                     col.ind = pcadf[,ncol(pcadf)],, # color by groups
                     palette = c("#999999", "#E69F00"), pointsize = 3,
                     addEllipses = TRUE, ellipse.type = "convex",
                     legend.title = "Protection-Status", title='') +theme(text = element_text(size = 10,family = 'Helvetica'),axis.text = element_text(size = 8,family = 'Helvetica'), legend.text = element_text(size = 20)
  )
  ggsave(filename = 'pca.svg',plot=manuscript_plot_pca,device = 'svg',path='../../data/plots/',dpi = 320)
  write.infile(df.pca, file.path( "../../data/pca_results/pca.csv"), sep = "\t")

}
col<-intersecting_colnames()
data<-read_data(timepoint,intersecting_colnames(500)) #for Day 7 and top 500 features
metadata<-fread('../../data/metadatafulltraining.csv')
data<-data %>% rownames_to_column(.,'AnimalID') %>% mutate(AnimalID = as.integer(AnimalID)) %>% left_join(dplyr::select(metadata,AnimalID),by='AnimalID',multiple = 'first')
data<-data %>% column_to_rownames(.,'AnimalID')
data$protectionStatus<- gsub('NonProt','non-protected',data$protectionStatus)
data$protectionStatus<- gsub('Prot','protected',data$protectionStatus)
pca_fun(data)
