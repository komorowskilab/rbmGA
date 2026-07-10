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

library(tidyverse)
library(data.table)
library(stringr)
library(R.ROSETTA)
library(csVisuNet)
library(RCy3)

read_data<-function(x,col){
  print(x)
  #reading the Day 7 data
  data<-fread(paste0('../../data/Rosetta_Decision_Tables/',x,'.csv'))

  #moving animal ids to rownames inorder to keep information intact
  data <- data %>% column_to_rownames(.,'AnimalID')

  #data preparation by subsetting the data with 70 genes and outcome
  df<-subset(data, select=c(col,'Group','protectionStatus'))
  df<-df[df$Group =='S' | df$Group =='O' | df$Group =='ABL1_a',]
  df<-df[df$protectionStatus=='Prot' | df$protectionStatus=='NonProt',]
  ros_data<- subset(df, select=c(col,'protectionStatus')) #

  return(ros_data)
}
intersecting_colnames<-function(val=500,timepoint='D7'){
  #MCFS selection
  res<-readRDS('../../data/Mcfs_results/mcfs_D7_gp1.rds')
  res_1<-readRDS('../../data/Mcfs_results/mcfs_D133_gp1.rds')
  sig_genes<- res$RI$attribute[1:val]
  sig_genes_1<- res_1$RI$attribute[1:val]
  union(sig_genes,sig_genes_1)
  intersect(sig_genes,sig_genes_1)
  col<-intersect(sig_genes,sig_genes_1)

  if(timepoint=='D0'){
    res<-readRDS('../../data/Mcfs_results/D0.rds')
    col<-res$RI$attribute[1:500]
  }

  return(col)
}

ros_data<-read_data('D7',intersecting_colnames(val=500))

rules<-readRDS('../../data/Rosetta_results/prediction_model.rds')


rules<- rules %>% distinct(features,levels,decision,.keep_all = T)

recal<-recalculateRules(ros_data,rules)

recal$decision<-gsub('NonProt','Non-protected',recal$decision,fixed = T)
recal$decision<-gsub('Prot','Protected',recal$decision,fixed = T)
recal$features<- gsub('X00000061294','IGHV3-43', recal$features,fixed=T)
recal$features<- gsub('LOC114669850','CD300LD', recal$features,fixed=T)
recal$features<- gsub('LOC702786','KRTAP16-1', recal$features,fixed=T)
recal$features<- gsub('KRTAP29.1','KRTAP29-1', recal$features,fixed=T)
recal$features<- gsub('X00000037843','OSCAR', recal$features,fixed=T)
recal$features<- gsub('C5H4orf45','SPMIP2', recal$features,fixed=T)
recal$features<- gsub('LOC701920','HMSD', recal$features,fixed=T)
recal$features<- gsub('LOC715623','AKR7A3', recal$features,fixed=T)
recal$features<- gsub('LOC715623','AKR7A3', recal$features,fixed=T)
recal$features<- gsub('X00000053012','LY9', recal$features,fixed=T)
recal$features<- gsub('LOC719667','SIGLEC6', recal$features,fixed=T)


cytoscapePing()

library(arcdiagram)

vis<-visunetcyto(recal,minAcc = 0.8,minSupp = 6,minDecisionCoverage = 0.1)

svg("../../data/plots/Cell Death and Inflammatory signalling.svg", width = 7, height = 7)
visuArc(vis,'Protected',feature = 'F2RL3',label='Cell Death and Inflammatory signalling')
dev.off()

svg("../../data/plots/Inflammation Control signalling.svg", width = 7, height = 7)
visuArc(vis,'Protected',feature = 'SPMIP2', label='Inflammation Control')
dev.off()





#View(recal[str_detect(recal$features,'C5H4orf45') & recal$decision=='Prot',] %>% distinct(features,levels,decision,.keep_all=TRUE))



