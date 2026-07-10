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
library(gridExtra)
#library(VisuNet)
library(import)
#import::from("~/Documents/gFunctions/weigthedpredictClass.R",predictClass)
import::from("../predictClass.R",predictClass)
import::from("../combi_gain.R",combi_gain)
#import::from("../visuArc.R",visuArc)

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

file_path='../../data/Rosetta_Decision_Tables/'
list_of_files <- list.files(path=file_path, pattern='.csv', all.files=FALSE, full.names=FALSE)
ros_data<-fread('../../data/Rosetta_Decision_Tables/D7.csv')
ros_data<-ros_data[ros_data$protectionStatus!='Uninfected',]
ros_data<- ros_data %>% column_to_rownames(.,'AnimalID')

classifier<-readRDS('../../data/Rosetta_results/prediction_model.rds')

col<-intersecting_colnames()
report<-data.frame(matrix(nrow=0,ncol=1))
names(report)<-'Accuracy'

groups<-c("ABL1_a" , "ABL2_A1" ,"ABL2_A2" ,"ABL2_A3" ,"ABL2_B",  "ABL3_A",  "ABL3_B",  "ABL3_C",  "O","S")
remove<- c('ABL1_a','O','S', "ABL3_B",  "ABL3_C","ABL2_A1" ,"ABL2_A2" ,"ABL2_A3")
groups<-groups[! groups %in% remove]

report_time<-data.frame()
report_dataset<-data.frame()

for(group in groups){
    predf<-ros_data[ros_data$Group %in% group,]
    predf<- subset(predf, select=c(col,'protectionStatus'))
    predf<-predf[predf$protectionStatus=='Prot' | predf$protectionStatus=='NonProt' ,]
    pred<-predictClass(predf[,1:(ncol(predf)-1)],classifier %>% distinct(features,levels,decision,.keep_all = T),,discrete = FALSE,validate = TRUE,defClass = predf$protectionStatus,normalizeMethod = 'mean',normalize = F)
    print(table(pred$out$predictedClass,pred$out$currentClass))
    print(pred)
    report[group,'Accuracy']<-pred$accuracy
    report_dataset[group,'Accuracy']<-pred$accuracy
  }



for (f in list_of_files){
  filename<-str_split(f,'[.]')[[1]][1]
  print(filename)
  predf<-fread(paste0('../../data/Rosetta_Decision_Tables/',f))
  predf<-predf[predf$Group =='S' | predf$Group =='O' | predf$Group =='ABL1_a',]
  predf<- predf %>% column_to_rownames(.,'AnimalID')
  predf<- subset(predf, select=c(col,'protectionStatus'))
  pred<-predictClass(predf[,1:(ncol(predf)-1)],classifier %>% distinct(features,levels,decision,.keep_all = T),,discrete = FALSE,validate = TRUE,defClass = predf$protectionStatus,normalize = F,weighted=TRUE)
  print(table(pred$out$predictedClass,pred$out$currentClass))
  print(pred)
  #print(pred$accuracy)

  if(filename=='D133'){
    assign('pred_D133',pred)
  }

  report[filename,'Accuracy']<-pred$accuracy
  report_time[filename,'Accuracy']<-pred$accuracy


}

report<- report %>% rownames_to_column(.,'Groups')
report$Groups<- factor(report$Groups, levels=c("ABL2_A1" ,"ABL2_A2" ,"ABL2_A3" ,"ABL2_B",  "ABL3_A",  "ABL3_B",  "ABL3_C","D0", "D3", "D7", 'D126','D129','D133','PreChal'))
#fwrite(report,'classifier_perf_on_groups.csv')

report_dataset<- report_dataset %>% rownames_to_column(.,'Groups')
report_dataset$Groups<- factor(report_dataset$Groups, levels=c("ABL2_A1" ,"ABL2_A2" ,"ABL2_A3" ,"ABL2_B",  "ABL3_A",  "ABL3_B",  "ABL3_C"))


report_time<- report_time %>% rownames_to_column(.,'Groups')
report_time$Groups<- factor(report_time$Groups, levels=c("D0", "D3", "D7", 'D126','D129','D133','PreChal'))



ggplot(report,aes(x=Groups,y=Accuracy))+geom_bar(stat="identity",position="dodge",fill='white',color='black')+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
xlab('External Datasets')+geom_text(aes(label = round(Accuracy,digits=2)),vjust=-2, color="black", size=5) +theme_bw() +
theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=1,size = 12))



x<-ggplot(report_time,aes(x=Groups,y=Accuracy))+geom_bar(stat="identity",position="dodge",fill='brown',color='black')+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  xlab('Partially External Datasets, Except D7 which is reclassification')+geom_text(aes(label = round(Accuracy,digits=2)),vjust=-2, color="black", size=5) +theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=1,size = 12))



y<-ggplot(report_dataset,aes(x=Groups,y=Accuracy))+geom_bar(stat="identity",position="dodge",fill='#69b3a2',color='black')+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
  xlab('External datasets; Day 7')+geom_text(aes(label = round(Accuracy,digits=2)),vjust=-2, color="black", size=5) +theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=1,size = 12))



grid.arrange(x,y,nrow=1)

net_classifier<-classifier

net_classifier$decision<-gsub('Prot','Protected',net_classifier$decision) #for paper network

# vis<-visunet(net_classifier[ , colSums(is.na(net_classifier)) == 0])
# visuArc(vis,decision = 'Protected',feature = 'F2RL3',label='Cell Death')

#fwrite(net_classifier[ , colSums(is.na(net_classifier)) == 0],'supp.csv',sep = '\t')



print(table(pred_D133$out$currentClass,pred_D133$out$predictedClass))

pred_D133$out$currentClass<-gsub('NonProt','non-protected',pred_D133$out$currentClass)
pred_D133$out$currentClass<-gsub('Prot','protected',pred_D133$out$currentClass)

pred_D133$out$predictedClass<-gsub('NonProt','non-protected',pred_D133$out$predictedClass)
pred_D133$out$predictedClass<-gsub('Prot','protected',pred_D133$out$predictedClass)

u <- union(pred_D133$out$predictedClass, pred_D133$out$currentClass)
t <- table(factor(pred_D133$out$predictedClass, u), factor(pred_D133$out$currentClass, u))



df<-caret::confusionMatrix(t)

as.data.frame(df$table) %>% ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  xlab('True Label') + ylab('Predicted label') +
  theme_bw() + theme(legend.position = "none") + coord_flip()

font<-ggplot2::theme_light(base_size = 12)

print(cvms::plot_confusion_matrix( as.data.frame(df$table),target_col = 'Var1', prediction_col = 'Var2', counts_col = 'Freq',   font_counts = cvms::font(
  size = 12,
  color = "black"
), add_row_percentages = F, add_col_percentages = T, add_normalized = F,font_col_percentages = cvms::font(size=5,color='black'))+ ggplot2::theme_light(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)))


