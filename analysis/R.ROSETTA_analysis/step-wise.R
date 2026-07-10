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
library(R.ROSETTA)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(randomForest)
library(caret)

#declare variables step-wise selection
step=10
timepoint<-'D7'
reducer <- 'Genetic'

#R.ROSETTA analysis on mcfs results
boost_check_case<- function(mcfs_result,col_size){
  if(mcfs_result$cutoff_value < 10){
    cat('Significant features from MCFS less than 10')
    cat('\nBoosting Features will be done in order of significance from MCFS')
    if(col_size > 500){
      case_num<-2
      cat('\nUpto 500 features')
    }else{case_num<-1
    cat('\nUpto number of features present in data')}
  } else if(mcfs_result$cutoff_value > 10 & mcfs_result$cutoff_value < 100 ){
    cat('Significant features from MCFS more than 10 but less than 100')
    cat('\nBoosting Features will be done in order of significance from MCFS')
    if(col_size > 500){
      case_num<-4
    }else{case_num<-3
    cat('\nUpto number of features present in data')}
  } else{
    cat('Significant features from MCFS more than 100')
    cat('\nBoosting Features will not be done')
    case_num<-5
  }
  return(case_num)
}
select_cols_boost<-function(boost_case_num,col_size,mcfs_result){
  if(boost_case_num==1 | boost_case_num==3){
    cols<- mcfs_result$RI[1:(col_size-1),'attribute']
  }else if(boost_case_num==2 | boost_case_num==4){
    cols<- mcfs_result$RI[1:500,'attribute']
  }else(boost_case_num==5)( cols<- mcfs_result$RI[1:mcfs_result$cutoff_value,'attribute'])
  return(cols)
}

mcfs_to_rosetta<-function(mcfs_df,out_name=NULL){

  response_var <- names(mcfs_df)[ncol(mcfs_df)]
  col_size<-ncol(mcfs_df)  #Should be in GLobal Env

  #requires boost_check_case function
  #requires select_cols_boost function
  #requires boost rosetta

  #Define paths and read MCFS file
  mcfs_file_path<- 'Mcfs_results/'

  if(is.null(out_name)){
    mcfs_result<-readRDS(paste0(mcfs_file_path,'/', list.files(path = 'Mcfs_results/',pattern = '.rds',recursive = FALSE))
    )
  }
  else{
    mcfs_file<- paste0(mcfs_file_path,out_name,'.rds')
    mcfs_result<-readRDS(mcfs_file)
  }


  #Identify boost case number
  boost_case_num<- boost_check_case(mcfs_result,col_size)
  cols<-select_cols_boost(boost_case_num,col_size,mcfs_result)

  #Make data for R:ROSETTA
  colnames(mcfs_df)<- gsub('-','_',colnames(mcfs_df))
  ros_data<- subset(mcfs_df,select=c(cols,response_var))
  #print(dim(ros_data))

  ros_results<-boost_ros(ros_data)
  #print(ros_results[3])
  cluster_df<-as.data.frame(ros_results[1])
  recal<-as.data.frame(ros_results[3])
  cluster_rules(cluster_df,recal)
}


# Use : boost_ros(ncol(autcon),autcon)
quality_check<- function(ros,quality){
  if(length(quality)==1){
    return (TRUE)
  }else if(ros$quality$accuracyMean > max(quality)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
boost_ros<-function(ros_data,step=2,discrete=FALSE,
                     reducer='Johnson',
                     JohnsonParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95)){
  library(ggplot2)
  ###requires quality check

  col_size<-ncol(ros_data)-1
  decision_var<-names(ros_data)[ncol(ros_data)]
  decision_roc<-names(table(eval(parse(text=paste("ros_data$", decision_var, sep = ""))))[1])
  df<-data.frame()

  print(col_size)
  count=10
  while(count<=plyr::round_any(col_size, step, f = ceiling) ){
    if(count>=col_size){
      print('this')
      ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,clroc='Prot',discrete = FALSE,discreteMethod = 'EqualFrequency',discreteParam = 3,discreteMask = FALSE,reducer = reducer,ruleFiltration = TRUE,cvNum = 10)
      quality<-append(quality,ros$quality$accuracyMean)

      ######################################
      if(quality_check(ros,quality)==TRUE){
        print('here1')
        ros_max<-ros
        data_max <-ros_data;
        recal<-recalculateRules(ros_data,ros$main,discrete = discrete)
      }
      ######################################

      response_var <- names(ros_data)[ncol(ros_data)]
      # Use '.' to represent all other columns as predictors
      formula <- as.formula(paste(response_var, "~ ."))

      #metrics
      df[col_size,'Nof']<-count
      df[col_size,'RBM_accuracy']<-ros$quality$accuracyMean
      df[col_size,'RBM_accuracy_std']<-ros$quality$ROC.AUC.MEAN

      #Counter
      count<-count+step
    }
    else{
      if(count==10)(quality<-c())
      trim_ros_data<-ros_data[,1:count]
      print(dim(trim_ros_data))
      trim_ros_data$decision<-eval(parse(text=paste("ros_data$", decision_var, sep = "")))
      trim_ros_data$decision<-factor(trim_ros_data$decision)
      ros<-rosetta(trim_ros_data,reducerDiscernibility = 'Object',roc=TRUE,clroc='Prot',discrete = FALSE,discreteMethod = 'EqualFrequency',discreteParam = 3,discreteMask = FALSE,reducer = reducer,ruleFiltration = TRUE,cvNum = 10)
      quality<-append(quality,ros$quality$accuracyMean)


      ######################################
      if(quality_check(ros,quality)==TRUE){
        print('here2')
        ros_max<-ros
        data_max <-trim_ros_data;
        recal<-recalculateRules(trim_ros_data,ros$main,discrete = discrete)
      }
      #####################################

      #metrics
      df[count,'Nof']<-count
      df[count,'RBM_accuracy']<-ros$quality$accuracyMean
      df[count,'RBM_accuracy_std']<-ros$quality$ROC.AUC.MEAN

      #Counter
      count<-count+step

    }
    print(quality)
  }
  return(list(data_max,ros_max,recal,df))
}
######################################## PROJECT DEPENDENT ################################################################

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

ros_data<-read_data(timepoint,intersecting_colnames(val=500))

#DISCRETIZATIOn

boost_res<-boost_ros(ros_data,step=step,reducer =reducer )

na.omit(as.data.frame(boost_res[4])) %>% ggplot()  +
   geom_point(aes(Nof,RBM_accuracy,shape='meanAccuracy',color='meanAccuracy'),size=5,,show.legend = T)+
   geom_point(aes(Nof,RBM_accuracy_std,shape='meanROC',color='meanROC'),size=5,show.legend = T) +
   scale_x_continuous(breaks = (as.data.frame(boost_res[4])$Nof))+   scale_y_continuous(breaks = seq(0,1,0.1)) + expand_limits(x = 10, y = 0) +
   scale_shape_manual(values = c("meanAccuracy" = 3, "meanROC" = 4)) +
   scale_color_manual(values = c("meanAccuracy" = "black", "meanROC" = "red")) +
   labs(x='Number of top features selected ', y='Accuracy', shape='Metric', color='Metric')+
   ggtitle('Step-wise Selection_Model') + #EDIT THIS To give title
   theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
   theme_minimal()

#
#   (as.data.frame(boost_res[4]))%>% ggplot(aes(Nof,RBM_accuracy)) + geom_line() + geom_point(size=4,color='black')+ggtitle('Boosting Performance') + geom_ribbon(aes(y = RBM_accuracy, ymin = RBM_accuracy - RBM_accuracy_std, ymax = RBM_accuracy + RBM_accuracy_std), alpha = .2) +
#     labs(x='Number of top features selected ', y='Accuracy')+ scale_x_continuous(breaks = (as.data.frame(boost_res[4])$Nof))+   scale_y_continuous(sec.axis = sec_axis(~ . / 10, name = "Conversion Rate (%)"))
#   theme_bw()  +
#     theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#
#
#
#   (as.data.frame(boost_res[4])) %>% ggplot(aes(Nof,RBM_accuracy))  + geom_point(size=5,color='black',shape=3,show.legend = T)+ geom_point(aes(Nof,RBM_accuracy_std),size=5,shape=4,color='red',show.legend = T) +
#     ggtitle('Step-wise Selection') + guides(color= guide_colourbar(position = "bottom")) +
#     labs(x='Number of top features selected ', y='Accuracy')+ scale_x_continuous(breaks = (as.data.frame(boost_res[4])$Nof))+   scale_y_continuous(breaks = seq(0,1,0.1)) + expand_limits(x = 10, y = 0) +
#     theme_bw()  +
#     scale_shape_manual(name = "Legend",
#                        labels = c("RBM_accuracy_std", "RBM_accuracy"),
#                        values = c(4, 3)) +
#     theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#
#
#
#
#
#
#
#
#
