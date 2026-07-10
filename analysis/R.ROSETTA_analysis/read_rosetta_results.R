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
#library(VisuNet)
library(import)
library(grid)
library(gridExtra)
#import::from('combi_gain.R',combi_gain)
#import::from('cluster_rules.R',cluster_rules)


#################################################### Read the rules from consensus model ####################################

rosetta_file_path<-'../../data/Rosetta_results/consensus_run/'
list_of_files<-list.files(rosetta_file_path,,pattern = 'ros_run*')
report_df<-data.frame()
nrls_man<-data.frame()
for (files in list_of_files){

  ros<-readRDS(paste0(rosetta_file_path,files))
  report_df[files,'Accuacy']<-ros$quality$accuracyMean
  nrls_man<- dplyr::bind_rows(nrls_man,ros$main)
}

consensus_rules<- nrls_man%>% distinct(features,levels,decision,.keep_all = T)


saveRDS(consensus_rules,'../../data/Rosetta_results/consensus_rules.rds') #save the rules from the model.



#################################################### Summarise the rules from model ####################################


table(consensus_rules$decision)
df<-as.data.frame(table(consensus_rules$decision,consensus_rules$accuracyRHS))
df$Var2<-as.numeric(as.character(df$Var2))
reshape2::melt(df %>%
  dplyr::mutate(condition_1 = Var2>=1,
                condition_2 = Var2>=0.9 & Var2<1,
                condition_3 = Var2<0.9) %>%
  group_by(Var1) %>%
  dplyr::summarise(
    'Accuracy=1 '= sum(Freq[condition_1] ),
    'Accuracy>0.9' = sum(Freq[condition_2] ),
    'Accuracy<0.9'= sum(Freq[condition_3])) ) %>% ggplot(aes(y=Var1, x=value,fill=variable)) + scale_fill_manual(name='Accuracy cutoffs',values = c("forestgreen", "yellow", "orange"))+ geom_bar(stat="identity",position="dodge") +xlab('No. of Rules') +ylab('Protection Status')+ theme_bw()
#%>% ggplot(aes(x=Var1,))





#################################################### PLOT ROC-CURVE FOR EACH CLASS ####################################


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


ros_data<-read_data('D7',intersecting_colnames(500))
ros_data<- ros_data %>% dplyr::select(1:50,protectionStatus) #since the best model was for 50 genes
roc_plot<-function(ROCstats,i){

  ROCstats<-ROCstats[ROCstats$CVNumber<10,]
  OMSpec <- rowMeans(unstack(ROCstats, form = OneMinusSpecificity ~ CVNumber))
  Sens <- rowMeans(unstack(ROCstats, form = Sensitivity ~ CVNumber))
  df<-data.frame('OMSpec'=OMSpec,'Sens'=Sens,'Class'=i)
  return(df)
}

#single run # ALSO plot ROC curve for each class
roc_df<-data.frame()
for ( i in  unique(ros_data[[ncol(ros_data)]])){
  ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,clroc=i,discrete = FALSE,discreteParam = 3,discreteMask = FALSE,reducer = 'Genetic',ruleFiltration = TRUE,cvNum = 10)
  roc_df<-roc_plot(ros$ROCstats,paste0(i,' - ', ros$quality$ROC.AUC.MEAN)) %>% bind_rows(roc_df)
}

if (!require(plotROC)) install.packages('plotROC')
library(plotROC)


#Ploting ROC curve
roc_df %>% ggplot(aes(x = OMSpec, y = Sens,colour = Class)) + scale_colour_manual(values=c("#999999", "#E69F00")) +
  geom_path() +
  geom_abline(lty = 5) +
  coord_equal() + style_roc(xlab = '1-specificity',ylab='Sensitivity')+
  theme_classic()+
  theme(legend.position = "bottom",legend.title = element_text(face = "bold"), legend.text = element_text(size = 8, colour = "black"))




