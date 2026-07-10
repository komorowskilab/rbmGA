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

#redundant for publication

intersecting_colnames<-function(val=500){
  #MCFS selection
  res<-readRDS('../data/Mcfs_results/mcfs_D7_gp1.rds')
  res_1<-readRDS('../data/Mcfs_results/mcfs_D133_gp1.rds')
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
  data<-fread(paste0('../data/Rosetta_Decision_Tables/',x,'.csv'))

  #moving animal ids to rownames inorder to keep information intact
  data <- data %>% column_to_rownames(.,'AnimalID')

  #data preparation by subsetting the data with 70 genes and outcome
  df<-subset(data, select=c(col,'Group','protectionStatus'))
  df<-df[df$Group =='S' | df$Group =='O' | df$Group =='ABL1_a',]
  df<-df[df$protectionStatus=='Prot' | df$protectionStatus=='NonProt',]
  ros_data<- subset(df, select=c(col,'protectionStatus')) #

  return(ros_data)
}

# cutoff_value<-function(val){
#   timepoints=c('D7','D133')
#
#   col<-intersecting_colnames(val=val)
#   #making the final decision table for R:ROSETTA
#   print(col)
#
#   for(x in timepoints){
#     #R:ROSETTA - NON linera modellung
#     ros_data<-read_data(x,col)
#     ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,clroc='Prot',discrete = FALSE,discreteMethod = 'EqualFrequency',discreteParam = 3,discreteMask = FALSE,reducer = 'Genetic',ruleFiltration = TRUE)
#     report[val,paste0('Mean_Accuracy_',x)]<-ros$quality$accuracyMean
#     report[val,paste0('Mean_ROC_Protected_',x)]<-ros$quality$ROC.AUC.MEAN
#     print(ros$quality)
#
#     ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,clroc='NonProt',discrete = FALSE,discreteMethod = 'EqualFrequency',discreteParam = 3,discreteMask = FALSE,reducer = 'Genetic',ruleFiltration = TRUE)
#     report[val,paste0('Mean_ROC_NonProtected_',x)]<-ros$quality$ROC.AUC.MEAN
#     print(ros$quality)
#   #  gc()
#   }
#
#  #
#  return(report)
#
# }
#
# report<-data.frame()
# for (i in seq(100,1000,50)){
#   report<-cutoff_value(i)
# }
#
# report<-na.omit(report)
# report<- report%>% rownames_to_column(.,'Number of features')
# fwrite(report,'../data/mcfs_report.csv')

report<-fread('../../data/mcfs_report.csv')
report<- report %>% column_to_rownames(.,'Number of features')
#plot report for Mean accuracy using R.ROSETTA
df_melt<- reshape2::melt(report,id.vars = report$`Number of features`)
df_melt$`Number of features`<-factor(df_melt$`Number of features`,levels=seq(100,1000,50))
df_melt[!str_detect(df_melt$variable,'ROC'),] %>% ggplot(aes(x=(`Number of features`),y=value)) +geom_line(aes(color = variable,group = variable)) +geom_point() + geom_vline(xintercept = as.factor(550), linetype="3313",  color = "blue", size=1.5)+ theme_bw()

#Day 7 and Day 133 results from MCFS
res<-readRDS('../../data/Mcfs_results/mcfs_D7_gp1.rds')
res_1<-readRDS('../../data/Mcfs_results/mcfs_D133_gp1.rds')

#RI_reportr for D7
RI_report<-data.frame()
for (i in seq(50,1000,50)){
  RI_report[i,'Sum']<-sum(res$RI[1:i,"RI"])
}
RI_report<-na.omit(RI_report)
RI_report <- RI_report %>% rownames_to_column(.,'Number of Features')
RI_report$`Number of Features`<- factor(RI_report$`Number of Features`,levels=seq(50,1000,50))
RI_report<- RI_report %>% mutate(pct_change = (Sum/lag(Sum) - 1))
RI_report[is.na(RI_report)]<-0

RI_report %>% ggplot(aes(x=`Number of Features`,y=Sum))+ geom_bar(stat="identity",position="dodge",fill='white',color='black',size=1) +geom_line(aes(x=`Number of Features`,y=pct_change,group=1),,stat="identity",color="red",size=1) + theme_bw()


#Interecting genes
report_gene<-data.frame()
for (i in seq(100,1000,50)){
  report_gene[i,'No. Intersecing genes']<-length(intersecting_colnames(i))
}
report_gene<-na.omit(report_gene)
report_gene <- report_gene %>% rownames_to_column(.,'Number of Features')
report_gene$`Number of Features`<- factor(report_gene$`Number of Features`,levels=seq(50,1000,50))
report_gene %>% ggplot(aes(x=`Number of Features`,y=`No. Intersecing genes`))+
geom_bar(stat="identity",position="dodge",fill='white',color='black',size=1)+
geom_text(aes(label=`No. Intersecing genes`), vjust=1.6, color="black", size=4.5,fontface='bold')+
theme_bw()+
theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



#RI plots

diff_RI_7<- c(0, abs(diff(res$RI$RI[seq(0,1000,50)], lag = 1)))
diff_RI_133 <- c(0, abs(diff(res_1$RI$RI[seq(0,1000,50)], lag = 1)))

# diff_RI_7<-res$RI$RI[seq(0,1000,50)]
# diff_RI_133<-res_1$RI$RI[seq(0,1000,50)]
diff_RI<-data.frame(RI_7=diff_RI_7,RI_133=diff_RI_133)
diff_RI$Number<-seq(50,1000,50)
colors<-c('RI_7'='red','RI_133'='blue')
p<-diff_RI[-1,] %>% ggplot(aes(x=Number))+ geom_line(aes(x=Number,y=RI_7,color='RI_7'),size=1)+ geom_line(aes(x=Number,y=RI_133,,color='RI_133'),size=1)+geom_point(aes(x=Number,y=RI_7,color='RI_7'),size=4) +
geom_point(aes(x=Number,y=RI_133,,color='RI_133'),size=4)+theme_bw()+scale_x_continuous(breaks=seq(50,1000,50)) +
xlab('Number of Features') + ylab('Difference in Relative importance %')+ guides(colour=guide_legend(title='RI'))+
scale_color_manual(values = colors)+
theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p+guides(fill = guide_legend(title = "LEFT"))


