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
library(csVisuNet)
library(import)
import::from('../change_gene_names.R',data_prep)
import::from('../predictClass.R',predictClass)

########################################################### CHANGE ##########################################################
#timepoint
TP<- 'D7'

#Metadata
metadata<-fread('../../data/external_validation/target_renamed_clinical_ortholog.csv')
#clincal ortho  ../data/target_renamed.csv
#O22 ../data/target_O22.csv
#O23 ../data/target_O23.csv


#normalised count matrix
count_data<- fread('../../data/external_validation/1.norm_matrix_clinicalortho_WT_revised.csv')
#clinical otho '../data/1.norm_matrix_clinicalortho_WT_revised.csv'
#O22 ../data/1.norm_matrix_O22_revised_20260127.csv
#O23 ../data/1.norm_matrix_all_O23_revised_20260127.csv


count_data<-data_prep(count_data,transpose = FALSE) #convertig enseml to symbol for RMs





#extra anlaysis can be ignored
#count_data<- fread('../data/raw/D7/1.norm_matrix.txt')
#names(count_data)[1]<- 'Gene.stable.ID'



#classifier
classifier <- readRDS('../../data/Rosetta_results/prediction_model.rds')


########################################################## ##########################################################


#Indexing from the metadata
decision_col_no <- which(str_detect(names(metadata), regex("outcome", ignore_case = TRUE)))
tp_col_no <- which(str_detect(names(metadata), regex("timepoint", ignore_case = TRUE)))
#sample_col_no <- which(str_detect(names(metadata), regex("sample", ignore_case = TRUE)))
sample_col_no <- which(grepl(rownames(count_data)[10],metadata))

#display distribution
print(table(metadata[[tp_col_no]],metadata[[decision_col_no]]))

#Identify animal outcomes
decisions <- unique(metadata[[decision_col_no]])

TP_metadata <- metadata[ metadata[[tp_col_no]]==TP & metadata[[decision_col_no]] %in% decisions & metadata$Sex=='m',]
#TP_metadata <- metadata[ metadata[[tp_col_no]]==TP & metadata[[decision_col_no]] %in% decisions & metadata$CMV.status=='positive' ,]
#TP_metadata <- metadata[ metadata[[tp_col_no]]==TP & metadata[[decision_col_no]] %in% decisions & metadata$CMV.status=='naive'  & metadata$Sex=='F',]
samples_2_trim<-TP_metadata[[sample_col_no]]

message('Done with metadata')


message('Starting count data')


TP_count_data <- count_data[rownames(count_data) %in% samples_2_trim,]
#TP_count_data <- TP_count_data %>% rownames_to_column(.,'Sample')  %>% dplyr::left_join(TP_metadata %>% dplyr::select(sample_col_no,decision_col_no),by=c('Sample')) %>% column_to_rownames(.,'Sample')

TP_count_data <- TP_count_data %>% rownames_to_column(.,'Sample')  %>% dplyr::left_join(TP_metadata %>% dplyr::select(sample_col_no,decision_col_no),by=c('Sample'='V1')) %>% column_to_rownames(.,'Sample')

message('Count data done, making predictions')

message('The Classifier labels are being changed from Protected to Non-protected to Controller and Progressor respectively')

# classifier$decision<-gsub('NonProt','Progressor',classifier$decision)
# classifier$decision<-gsub('Prot','Controller',classifier$decision)

classifier$decision<-gsub('NonProt','nonprotected',classifier$decision)
classifier$decision<-gsub('Prot','protected',classifier$decision)

classifier$features<-gsub('C5H4orf45','SPMIP2',classifier$features)

i<-ncol(TP_count_data)
pred<- predictClass(TP_count_data,classifier %>% distinct(features,levels,decision,.keep_all = T),normalize = F,discrete = F,weighted = T, validate= T, defClass = TP_count_data[[i]])

message(paste0('Accuracy = ',pred$accuracy))

print(table(pred$out$currentClass,pred$out$predictedClass))

pred$out$currentClass<-gsub('nonprotected','non-protected',pred$out$currentClass)
pred$out$predictedClass<-gsub('nonprotected','non-protected',pred$out$predictedClass)

u <- union(pred$out$predictedClass, pred$out$currentClass)
t <- table(factor(pred$out$predictedClass, u), factor(pred$out$currentClass, u))



df<-caret::confusionMatrix(t)

as.data.frame(df$table) %>% ggplot(mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  xlab('True Label') + ylab('Predicted label') +
  theme_bw() + theme(legend.position = "none") + coord_flip()


cvms::plot_confusion_matrix( as.data.frame(df$table),target_col = 'Var1', prediction_col = 'Var2', counts_col = 'Freq',   font_counts = cvms::font(
  size = 12,
  color = "black"
), add_row_percentages = F, add_col_percentages = T, add_normalized = F,font_col_percentages = cvms::font(size=5,color='black'))+ ggplot2::theme_light(base_size = 20) + theme(plot.title = element_text(hjust = 0.5))

