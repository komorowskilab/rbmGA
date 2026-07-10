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
set.seed(1234)
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
library(csVisuNet)
library(RCy3)


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
  df<-df[df$protectionStatus=='Prot' | df$protectionStatus=='NonProt',]
  ros_data<- subset(df, select=c(col,'protectionStatus')) #

  return(ros_data)
}


#sorted genes into different files
cell_sort_path<-'../../data/Cell_type/'
#list_of_cell_files<-list.files(cell_sort_path,pattern = '.*phil')
#list_of_cell_files<- append(list_of_cell_files,list.files(cell_sort_path,pattern = 'monocyte'))
list_of_cell_files<-list.files(cell_sort_path,pattern = '.csv')

#removing them because they are not relevant or for instance of T-cell they do not form any rules
remove<- c('hepatocytes.csv','neutro_mono.csv','tcell.csv')
list_of_cell_files<-list_of_cell_files[! list_of_cell_files %in% remove]


file_path='files/'
list_of_files <- list.files(path=file_path, pattern='.csv', all.files=FALSE, full.names=FALSE)

report<-data.frame(matrix(nrow=0,ncol=1))
names(report)<-'Accuracy'

###UNCOMMENT TO GENERATE THEM AGAIN
fin_recal<-data.frame()
  for (cell_type in list_of_cell_files){
    print(cell_type)
    ros_data<-read_data('D7',intersecting_colnames(500))
    file_name<-str_split(cell_type,'[.]')[[1]][1]
    cs<-fread(paste0('~/Documents/Gale_lab/Data/Cell_type/',cell_type))
    cols<-cs$Current_SYMBOL
    # filename<-str_split(,'[.]')[[1]][1] #deveopment code
    # ros_data<-fread(paste0('files/',f)) #deveopment code
    ros_data<- subset(ros_data, select=c(cols,'protectionStatus'))
    ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,clroc='Prot',discrete = FALSE,discreteMethod = 'EqualFrequency',discreteParam = 3,discreteMask = FALSE,reducer = 'Genetic',ruleFiltration = TRUE)
    recal<-recalculateRules(ros_data,ros$main)
    recal<-recal[recal$supportLHS>6,]
    #fin_recal<-rbind.fill(fin_recal,recal)
    fin_recal<-dplyr::bind_rows(fin_recal,recal)
    report[file_name,'Accuracy']<-ros$quality$accuracyMean
  }


fin_recal<-fin_recal %>% group_by(features, levels, decision) %>% dplyr::slice(1, 2)%>%ungroup()


# report<- report %>% rownames_to_column(.,'Cell Types')
# a<-ggplot(report,aes(x=`Cell Types`,y=Accuracy))+geom_bar(stat="identity",position="dodge",fill='white',color='black')+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
#   xlab('Cross-validation scores for Cell Types')+geom_text(aes(label = round(Accuracy,digits=2)),vjust=-2, color="black", size=5) +theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=1,size = 8))
#
# print(a)
#save rules
#saveRDS(fin_recal,'../../data/Rosetta_results/immune_cell_type_classifier.rds')
#fin_recal<-readRDS('../../data/Rosetta_results/immune_cell_type_classifier.rds')
#cytoscapePing()

#CHANGING GENE NAMES
fin_recal$decision<-gsub('NonProt','Non-protected',fin_recal$decision,fixed = T)
fin_recal$decision<-gsub('Prot','Protected',fin_recal$decision,fixed = T)
fin_recal$features<- gsub('X00000061294','IGHV3-43', fin_recal$features,fixed=T)
fin_recal$features<- gsub('LOC114669850','CD300LD', fin_recal$features,fixed=T)
fin_recal$features<- gsub('LOC702786','KRTAP16-1', fin_recal$features,fixed=T)
fin_recal$features<- gsub('KRTAP29.1','KRTAP29-1', fin_recal$features,fixed=T)
fin_recal$features<- gsub('X00000037843','OSCAR', fin_recal$features,fixed=T)
fin_recal$features<- gsub('C5H4orf45','SPMIP2', fin_recal$features,fixed=T)
fin_recal$features<- gsub('LOC701920','HMSD', fin_recal$features,fixed=T)
fin_recal$features<- gsub('LOC715623','AKR7A3', fin_recal$features,fixed=T)
fin_recal$features<- gsub('LOC715623','AKR7A3', fin_recal$features,fixed=T)
fin_recal$features<- gsub('X00000053012','LY9', fin_recal$features,fixed=T)
fin_recal$features<- gsub('LOC719667','SIGLEC6', fin_recal$features,fixed=T)


#visunetcyto(fin_recal,minAcc = 0,minSupp =0,minDecisionCoverage = 0) #network visualisation

#check unique rules from main classifier
ros<-readRDS('../../data/Rosetta_results/consensus_rules.rds')
unique_rules<-dplyr::bind_rows(fin_recal,ros)
unique_rules<-distinct(unique_rules)
