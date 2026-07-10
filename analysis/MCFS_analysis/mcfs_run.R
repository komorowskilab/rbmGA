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
#load librarires
if (!require(data.table)) install.packages('data.table')
library(data.table)
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)
if (!require(rmcfs)) install.packages('rmcfs')
library(rmcfs)

#MCFS function and their dependencies
#function to check if the filename exist then create a new one to avoid overalpping
check_out_file<-function(){
  if(file.exists("Mcfs_results/out.txt")){
    cat('Creating new out file since out.txt already exists')
    new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'/Mcfs_results/'))
    file.create(paste0(getwd(),'/Mcfs_results/',new_indexed_file))
    return(paste0('Mcfs_results/',new_indexed_file))
  }
  else{
    new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'/Mcfs_results/'))
    file.create(paste0(getwd(),'/Mcfs_results/',new_indexed_file))
    return(paste0('Mcfs_results/',new_indexed_file))
  }
}
#check if the file is exist then create a new index for the file
create_indexed_file = function(filename,extension='.txt', file_path){

  if(!file.exists(paste0(file_path,filename,extension))){return(paste0(filename,extension))}
  i=1
  repeat {
    f = paste0(filename,i,extension)
    if(!file.exists(paste0(file_path,f))){return(f)}
    i=i+1
  }
}
#mcfs run
mcfs_run<-function(mcfs_df,out_name,permutations=3,projections=projections, projectionSize=projectionSize,splitSetSize=splitSetSize ,set_response_var=FALSE,response_var=NULL,plot_res=TRUE){
  library(rJava)
  library(rmcfs)

  # requires create_indexed_file
  # requires check_out_file
  ## mcfs_out_path : file path for saving mcfs file
  ## plot_out_path : file path for saving mcfs plots and graphs

  cat("Checking input arguments...\n")
  if(is_empty(out_name)){
    mystop('output name is empty')
  }
  cat('Checking for out file')
  out_file<-check_out_file()

  mcfs_df<-fix.data(mcfs_df)

  #define formula for mcfs call
  if(set_response_var==TRUE){
    response_var<-response_var
  }
  else{
    #set the response variable
    cat('Choosing the last column as Response variable')
    response_var <- names(mcfs_df)[ncol(mcfs_df)]
    # Use '.' to represent all other columns as predictors
    formula <- as.formula(paste(response_var, "~ ."))
  }

  #Response variable
  cat(paste0('The Response varibale is : ',response_var))
  cat(paste0('\nStarted running Monte Carlo Feature Selection at ',Sys.time()))

  #run monte carlo feature
  mcfs_result <- mcfs(formula, mcfs_df, projections='auto', projectionSize='auto', splits=5, splitSetSize = 500, cutoffPermutations=permutations,threadsNumber= 8,balance=2)

  #Save mcfs in R data format
  #create directory
  dir.create(file.path(getwd(), 'Mcfs_results'),showWarnings = FALSE)
  mcfs_out_path<- file.path(getwd(), 'Mcfs_results')
  saveRDS(mcfs_result,paste0(mcfs_out_path,'/',out_name,'.rds'))

  #Plot Graphs and results
  if(plot_res==TRUE){
    dir.create(file.path(getwd(), 'Mcfs_results/plots'),showWarnings = FALSE)
    plot_out_path<-file.path(getwd(), 'Mcfs_results/plots')
    gid <- build.idgraph(mcfs_result, size = 6, size_ID = 12, orphan_nodes = TRUE)
    export.plots(mcfs_result, label=out_name,mcfs_df, idgraph = gid, path = plot_out_path, color = "darkgreen",plot_format = 'png',size = 20,image_width = 15,image_height = 15,cex=2)
    dev.off() #Not sure if this is the right place but it works
  }

  #Information for reproducibility
  cat(paste0(out_name,':Mcfs file Saved in Directory : ',paste0(mcfs_out_path,'/',out_name,'.rds')))
  cat(paste0('\nMCFS run exited with response variable: ',response_var))
  cat('\nCheck Out file in directory for execution IFNO')
  cat('######## Monte Carlo Feature Selection ############',file=out_file,sep="\n",append=TRUE)
  cat(paste0('Files and plots saved in Directory : ',mcfs_out_path),file=out_file,sep="\n",append=TRUE)
  cat(capture.output(print.mcfs(mcfs_result)),file = out_file,sep="\n",append=TRUE)
}

#data files
file_path='../../data/Rosetta_Decision_Tables/'
list_of_files <- list.files(path=file_path, pattern='.csv', all.files=FALSE, full.names=FALSE)

#loop to run mcfs for all timepoints
for (f in list_of_files){
  filename<-str_split(f,'[.]')[[1]][1] #creating filename for mcfs outfiles
  data<-fread(paste0(file_path,f),encoding = 'UTF-8')
  data <- data %>% column_to_rownames(.,'AnimalID') #movind animal ids to rownames
  #set the experminetal parameters orelse it was set to default values
  mcfs_run(data,filename,permutations = 20,projections = 20000,projectionSize = 100, splitSetSize = 85)
  }
