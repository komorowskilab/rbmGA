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
library(org.Hs.eg.db)
library(R.ROSETTA)
library(MVN)
library(ggplot2)
library(stringr)
library(pheatmap)
library(funModeling)
#library(VisuNet)
library(plotROC)


#MAke directories ; data and src.
#store the code in src
#in folder data make a new directory named Rosetta_results it is case sensitive
#then run the program to achieve no error execution

#in future we can automate the directory process


create_indexed_file = function(filename,extension='.txt', file_path){

  if(!file.exists(paste0(file_path,filename,extension))){return(paste0(filename,extension))}
  i=1
  repeat {
    f = paste0(filename,i,extension)
    if(!file.exists(paste0(file_path,f))){return(f)}
    i=i+1
  }
}


check_out_file<-function(){
  if(file.exists("../../data/Rosetta_results//out.txt")){
    cat('Creating new out file since out.txt already exists')
    new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'../../data/Rosetta_results/'))
    file.create(paste0('../../data/Rosetta_results/',new_indexed_file))
    return(paste0('../../data/Rosetta_results/',new_indexed_file))
  }
  else{
    new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'../../data/Rosetta_results/'))
    file.create(paste0('../../data/Rosetta_results/',new_indexed_file))
    return(paste0('../../data/Rosetta_results/',new_indexed_file))
  }
}


genetic_consensus<-function(ros_data,num_of_iter=100,discrete=FALSE){
  x<-new.env(hash = TRUE, parent = parent.frame(), size = 29L)
  decision_var<-names(ros_data)[ncol(ros_data)]
  decision_roc<-names(table(eval(parse(text=paste("ros_data$", decision_var, sep = ""))))[1])
  cat('Creating Directory for consensus. results')
  dir.create(file.path('../../data/Rosetta_results/consensus_run'),showWarnings = FALSE)
  cat('Checking for out file')
  out_file<-check_out_file()
  cat('you can find the list of accepted runs in a RData file')
  accepted_list<-vector()

  for(i in 1:num_of_iter){
    set.seed(i) #to generate random initialization and for debugging
    #shuffle column here
    start_time<- Sys.time()
    print(start_time)
    print(paste('Iteration Number:',i))
    ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,fallBackClass = 'Prot',discrete = discrete,discreteMethod = 'EqualFrequency',reducer = 'Genetic',clroc=decision_roc,ruleFiltration = TRUE,GeneticParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple"))
    # put it into environment and then remove when merged to save memory #suggestion

    if(ros$quality$accuracyMean>0.790){
      assign(paste0('ros_',i),ros$main )
      accepted_list<- append(accepted_list,i)
    }else(cat(paste0('Rosetta result from run ',i,' is discarded from merging'),file = out_file,append=TRUE,sep="\n"))

    saveRDS(ros,paste0('../../data/Rosetta_results/consensus_run/ros_run_',i,'.rds'))

    end_time<- Sys.time()
    cat(paste('This iteration took in seconds',end_time-start_time),file = out_file,append=TRUE,sep="\n")

    if(i == 1){
      total_time<- ((end_time-start_time)/60)*num_of_iter
      cat(paste('Approximate time for execution of Consensus Genetic would be',total_time,'Minutes\n'),file = out_file,append=TRUE,sep="\n")
    }

  }

  saveRDS(accepted_list,'../../data/Rosetta_results/consensus_run/accepted_list.rds')
  if(!is_empty(accepted_list)){
    rbms<- mget(ls(pattern = "ros_[0-9]"))

    if(length(accepted_list)>1){
      nrls<-mergeRBMs(rbms,defClass = eval(parse(text=paste("ros_data$", decision_var, sep = ""))))
      return(nrls)
    }else(return(ros))

  }else(return(ros))

}


intersecting_colnames<-function(){
  #MCFS selection
  res<-readRDS('../../data/Mcfs_results/mcfs_D7_gp1.rds')
  res_1<-readRDS('../../data/Mcfs_results/mcfs_D133_gp1.rds')
  sig_genes<- res$RI$attribute[1:500]
  sig_genes_1<- res_1$RI$attribute[1:500]
  union(sig_genes,sig_genes_1)
  intersect(sig_genes,sig_genes_1)
  col<-intersect(sig_genes,sig_genes_1)
  return(col)
}




######################################## PROJECT DEPENDENT ################################################################

#reading the Day 7 data #YOU WILL CHANGE THE DF HERE
data<-fread('../../data/Rosetta_Decision_Tables/D7.csv')

#moving animal ids to rownames inorder to keep information intact
data <- data %>% column_to_rownames(.,'AnimalID')

#selecting the 70 genes from MCFS
col<-intersecting_colnames()

#data preparation by subsetting the data with 70 genes and outcome
df<-subset(data, select=c(col,'Group','protectionStatus'))
df<-df[df$Group =='S' | df$Group =='O' | df$Group =='ABL1_a',]
df<-df[df$protectionStatus=='Prot' | df$protectionStatus=='NonProt',]


#No of features

#making the final decision table for R:ROSETTA
ros_data<- subset(df, select=c(col[1:50],'protectionStatus')) #50 genes because of step-wise predictor selection

######################################## ######################################## ########################################

##IMPORTANT####
#running consensus genetic model
nrls<-genetic_consensus(ros_data = ros_data,num_of_iter = 600)
saveRDS(nrls,'../../data/Rosetta_results/nrls_1.rds')








