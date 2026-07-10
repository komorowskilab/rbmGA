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
#library(org.Hs.eg.db)
library(R.ROSETTA)
library(MVN)
library(ggplot2)
library(stringr)
library(pheatmap)
library(funModeling)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(RColorBrewer)
import::from("../cluster_rules.R",cluster_rules)
import::from("../rule_characterization.R",rule_characterization)

#timepoint
TP<-'D7'

read_data<-function(x,col,group=NULL){
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


groups<-list(c('ABL1_a','O','S'))
for (i in 1:length(groups)){
  print(groups[i])
  ros_data<-read_data(TP,intersecting_colnames(500),group=groups[i])
  classifier<-readRDS('../../data/Rosetta_results/prediction_model.rds')
  classifier<- classifier %>% dplyr::distinct(features,levels,decision,.keep_all = T)
  #classifier<-classifier[str_count(classifier$features,',')>0,]
  recal<-recalculateRules(ros_data,classifier)
  recal$decision<-as.character(recal$decision)
  mod_recal<-recal[(str_detect(recal$features,'F2RL3') | str_detect(recal$features,'C5H4orf45')) | str_detect(recal$features,'KRTAP29.1') & recal$decision=='Prot',]
  #mod_recal<-mod_recal %>% mutate(Cluster = case_when(str_detect(mod_recal$features,'F2RL3')~ "Cell_Death",str_detect(mod_recal$features,'C5H4orf45')~ "Inf_Control"),str_detect(mod_recal$features,'KRTAP29.1') & str_detect(mod_recal$features,'LOC702786')~ "Keratine" )

  mod_recal<-recal %>% mutate(Rules= case_when(str_detect(recal$features,'F2RL3')~ "Cell_Death",str_detect(recal$features,'C5H4orf45')~ "Inf_Control",str_detect(recal$decision,'NonProt')~ "NonProt",.default = 'Other_in_Prot'))
  assign(paste0('val',i),rule_characterization(ros_data,mod_recal,group=groups[i],type='Rules',TP=TP))



}

grid.arrange(val1[[1]],val1[[2]],nrow=1)








