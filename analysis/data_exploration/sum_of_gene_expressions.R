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
library(import)
library(gridExtra)
library(ggpubr)
library(ggrepel)


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

#leukocytereceptor cluster
LRC<-c('VSTM1','TARM1','LILRB1','A1BG','KIR3DL0')

#D7
df<-read_data('D7',intersecting_colnames(500))
df$Sum_of_expressions<-rowSums(df[,1:70])
p1<-df %>% ggplot(aes(x=protectionStatus,y=Sum_of_expressions,fill=protectionStatus)) +  geom_boxplot()+  geom_jitter(shape=16, position=position_jitter(0.2),size=3) +
  stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) +
scale_fill_manual(values=c("#999999", "#E69F00"))+ scale_y_continuous(limits = c(0,80),breaks=seq(0,80,20))+ theme_bw() + ggtitle(('Day 7')) + ylab('Sum of 70 gene expression values')+xlab('Protection Status')+
theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

t.test(df[df$protectionStatus=='Prot',"Sum_of_expressions"],df[df$protectionStatus=='NonProt',"Sum_of_expressions"],var.equal = F)

print(p1) #day 7




#D133
df<-read_data('D133',intersecting_colnames(500))
df$Sum_of_expressions<-rowSums(df[,1:70])
p2<-df %>% ggplot(aes(x=protectionStatus,y=Sum_of_expressions,fill=protectionStatus)) +  geom_boxplot()+  geom_jitter(shape=16, position=position_jitter(0.2),size=3) +
  stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) +
scale_fill_manual(values=c("#999999", "#E69F00"))+ scale_y_continuous(limits = c(0,80),breaks=seq(0,80,20))+ theme_bw() + ggtitle(('Day 133')) + ylab('Sum of 70 gene expression values')+xlab('Protection Status')+
theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

t.test(df[df$protectionStatus=='Prot',"Sum_of_expressions"],df[df$protectionStatus=='NonProt',"Sum_of_expressions"],var.equal = F)

print(p2) #day 133
