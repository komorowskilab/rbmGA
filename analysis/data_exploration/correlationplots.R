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
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(plyr)
library(rmcfs)
library(edgeR)
library(sva)
library(parallel)
library(doParallel)
library(R.ROSETTA)
library(VisuNet)
library(import)
library(corrplot)
library(RColorBrewer)
library(MVN)
library(gridGraphics)
library(gridExtra)

timepoint<-'D7' #you can change it to 'D133'

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


#read data
read_data<-function(x,col,decision=NULL){
  print(x)
  #reading the Day 7 data
  data<-fread(paste0('../../data/Rosetta_Decision_Tables/',x,'.csv'))

  #moving animal ids to rownames inorder to keep information intact
  data <- data %>% column_to_rownames(.,'AnimalID')

  #data preparation by subsetting the data with 70 genes and outcome
  df<-subset(data, select=c(col,'Group','protectionStatus'))
  df<-df[df$Group =='S' | df$Group =='O' | df$Group =='ABL1_a',]
  df<-df[df$protectionStatus==decision,]
  ros_data<- subset(df, select=c(col,'protectionStatus')) #

  return(ros_data)
}



analysis<-function(data,decision){
  #Multivariate analysis
  #result <- MVN::mvn(data = data[,1:(ncol(data)-1)], mvnTest = "royston",univariateTest = "SW", desc = TRUE, univariatePlot = 'qqplot', showOutliers = TRUE,showNewData = TRUE)
  result<-MVN::mvn(data = data[,1:70], mvnTest = "royston", desc = TRUE, showOutliers = TRUE,showNewData = TRUE, multivariatePlot = 'qq') #1:70 because lst column is their decision
  grid.echo()
  result <- grid.grab()
  #plot(result, type = "qq")

  #correlation analysis
  M <-cor(data[,1:70])
  scalebluered <- colorRampPalette(brewer.pal(8, "RdBu"))(8)
  corr<-corrplot(M, order="hclust",col=scalebluered, tl.cex = 0.7,tl.col = 'grey39',tl.srt=45,type='lower',method = 'ellipse',title = paste0(decision,'ected'),mar=c(0,0,2,0))
  grid.echo()
  corr <- grid.grab()
  return(list(result,corr))
}

decision<-c('Prot','NonProt')
for(i in decision){
  data<-read_data(timepoint,intersecting_colnames(500),decision=i) #for Day 7 and top 500 features
  res<-analysis(data,i)
  assign(paste0('x_',i),res)
}

#grid.arrange(x_Prot[[2]],x_NonProt[[2]],nrow=1) # if you want to arrange them, might overlap
