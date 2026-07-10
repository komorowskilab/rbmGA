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

library(R.ROSETTA)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(randomForest)
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

shuffle_data<-function(ros_data){
  df_shuffled=transform( ros_data, protectionStatus = sample(protectionStatus))
  return(df_shuffled)
}


######################################## PROJECT DEPENDENT ################################################################
ros_data<-read_data('D7',intersecting_colnames(500))


######################################## ######################################## ########################################


sink(file = "lm_output.txt")
print('Orginal labels')
print(ros_data$protectionStatus)

print('Now the shuffling starts')
quaility<-data.frame()

#THINK BEFORE RUNNING
# for(i in 1:1000){
#   set.seed(i) #to generate random initialization and for debugging
#   df<-shuffle_data(ros_data)
#   df<- df %>% dplyr::select(1:50,'protectionStatus') #because 50 genes model generated the best model
#   print(df$protectionStatus)
#   ros<-rosetta(df,reducerDiscernibility = 'Object',roc=TRUE,clroc='Prot',discrete = FALSE,discreteMethod = 'EqualFrequency',discreteParam = 3,discreteMask = FALSE,reducer = 'Genetic',ruleFiltration = TRUE,cvNum = 10)
#   print(ros$quality$accuracyMean)
#   quaility[i,'Accuracy']<-ros$quality$accuracyMean
#
# }
sink(file = NULL)
#saveRDS(quaility,'quality.rds')

#directly reading the quality values produced. You can run the code above and comment the line below to reproduce the results.
quaility<-readRDS('../../data/Rosetta_results/quality.rds')

library(grid)
grob <- grobTree(textGrob("n=1000", x=0.8,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=13, fontface="italic")))



ci_bounds <- quantile(quaility$Accuracy, probs = c(0.025, 0.975), na.rm = TRUE)

print(quaility %>%
  ggplot(aes(x = Accuracy)) +
  geom_histogram(aes(y = ..density..), fill = 'royalblue1', color = "white", alpha = 0.8) +
  geom_density(color = "black", size = 1) +
  geom_vline(xintercept = ci_bounds, color = "firebrick", linetype = "dotted", size = 1) +
  geom_vline(aes(xintercept = 0.9, color = "p<0.001"), size = 1.5, linetype = "longdash") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  annotation_custom(grob) +
  guides(colour = guide_legend(title = 'Original Classifer'))+
  annotate("text", x = 0.70, y = 1, label = "Distribution of \nRandom Classifiers", fontface = "italic", color = "darkblue", size = 4) +
  annotate("segment", x = 0.7, y = 0.8, xend = 0.65, yend = 0.4,  arrow = arrow(length = unit(0.2, "cm")), color = "darkblue"))



# quaility %>% ggplot() + geom_histogram(aes(x=Accuracy),fill='royalblue1') + scale_x_continuous(breaks = seq(0,1,0.1))+geom_vline(aes(xintercept=0.9,color="Orginal Classifier: p<0.001"),size = 1.5, linetype = "longdash") +
#   theme_classic() + annotation_custom(grob)+ guides(colour=guide_legend(title='Original Classif'))






