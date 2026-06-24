#' @title Consensus Run
#' @author Girish Pulinkala
#' @description multiple runs of RBM to reduce stocashicity of genetic reducer
#' @import tidyverse
#' @import data.table
#' @import tibble
#' @import stringr
#' @import dplyr
#' @import R.ROSETTA
#'
#' @export
#'
#' @param data decision table used in the training
#' @param num_of_iter number of iterations
#' @param reducer A character containing name of reducer method: Johnson or Genetic. Default is Johnson.
#' @param discrete Logical. Set TRUE for discrete data. Default is FALSE.
#' @param discreteMethod A character containing discretization method: EqualFrequency, MDL, Naive, SemiNaive or BROrthogonal. Default is EqualFrequency.
#' @param discreteParam A vector containing discretization parameters. May be of different length and values. See examples.
#' @param discreteMask Logical. Set FALSE to disable discretization mask. Default is TRUE.
#' @param ruleFiltration Logical. Set TRUE to filter out rules. Default is TRUE. UNLIKE Original R.ROSETTA code
#' @param JohnsonParam A vector containing Johnson reducer parameters.
#' @param GeneticParam A vector containing Genetic reducer parameters. description
#'
#' @examples
#' # example code
#' #rules<-consensus_run(autcon)
#'




consensus_run<-function(data,num_of_iter=100,reducer = 'Johnson',discrete=FALSE,discreteMethod = 'EqualFrequency',
                        cvNum=10,ruleFiltration = TRUE,discreteParam =3,
                        JohnsonParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95),
                        GeneticParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple") ){

  #foundational functions
  make_directory <- function(file_path){

    message('Creating Directory for consensus results')
    message(paste0('this is you will find your results:', file_path,'/data/Rosetta_results/consensus_run'))
    dir.create(file.path(paste0(file_path,'/data/Rosetta_results/consensus_run')),showWarnings = TRUE,recursive = TRUE)
    message('Checking for out file')
    out_file<-check_out_file(file_path)
    message('you can find the list of accepted runs in a RData file')
    return(out_file)
  }


  create_indexed_file = function(filename,extension='.txt', file_path){

    if(!file.exists(paste0(file_path,filename,extension))){return(paste0(filename,extension))}
    i=1
    repeat {
      f = paste0(filename,i,extension)
      if(!file.exists(paste0(file_path,f))){return(f)}
      i=i+1
    }
  }


  check_out_file<-function(file_path){
    if(file.exists(paste0(file_path,"/data/Rosetta_results//out.txt"))){
      cat('Creating new out file since out.txt already exists')
      new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'/data/Rosetta_results/'))
      file.create(paste0(file_path,'/data/Rosetta_results/',new_indexed_file))
      return(paste0(file_path,'/data/Rosetta_results/',new_indexed_file))
    }
    else{
      new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'/data/Rosetta_results/'))
      file.create(paste0(file_path,'/data/Rosetta_results/',new_indexed_file))
      return(paste0(file_path,'/data/Rosetta_results/',new_indexed_file))
    }
  }

  #Main script starts

  x<-new.env(hash = TRUE, parent = parent.frame(), size = 29L)
  decision_var<-names(data)[ncol(data)]
  decision_roc<-names(table(eval(parse(text=paste("data$", decision_var, sep = ""))))[1])
  file_path<-getwd()
  out_file<-make_directory(file_path)
  accepted_list<-vector()

  for(i in 1:num_of_iter){
    set.seed(i) #to generate random initialization and for debugging
    #shuffle column here
    start_time<- Sys.time()
    print(start_time)
    print(paste('Iteration Number:',i))
    ros<-rosetta(data,reducerDiscernibility = 'Object',fallBackClass = decision_roc,discrete = discrete,discreteMethod = discreteMethod,discreteParam =discreteParam ,reducer = reducer,clroc=decision_roc,ruleFiltration = ruleFiltration,GeneticParam = GeneticParam, JohnsonParam = JohnsonParam,cvNum = cvNum )
    # put it into environment and then remove when merged to save memory #suggestion

    if(ros$quality$accuracyMean>0.790){
      assign(paste0('ros_',i),ros$main )
      accepted_list<- append(accepted_list,i)
    }else(cat(paste0('Rosetta result from run ',i,' is discarded from merging'),file = out_file,append=TRUE,sep="\n"))

    saveRDS(ros,paste0(file_path,'/data/Rosetta_results/consensus_run/ros_run_',i,'.rds'))

    end_time<- Sys.time()
    cat(paste('This iteration took in seconds',end_time-start_time),file = out_file,append=TRUE,sep="\n")

    if(i == 1){
      total_time<- ((end_time-start_time)/60)*num_of_iter
      cat(paste('Approximate time for execution of Consensus Genetic would be',total_time,'Minutes\n'),file = out_file,append=TRUE,sep="\n")
    }

  }

  saveRDS(accepted_list,paste0(file_path,'/data/Rosetta_results/consensus_run/accepted_list.rds'))

  if(!is_empty(accepted_list)){
    rbms<- mget(ls(pattern = "ros_[0-9]"))

    if(length(accepted_list)>1){

      rosetta_file_path<-paste0(file_path,'/data/Rosetta_results/consensus_run/')
      list_of_files<-list.files(rosetta_file_path,,pattern = 'ros_run*')
      report_df<-data.frame()
      nrls_man<-data.frame()
      for (files in list_of_files){

        ros<-readRDS(paste0(rosetta_file_path,files))
        #report_df[files,'Accuracy']<-ros$quality$accuracyMean
        nrls_man<- dplyr::bind_rows(nrls_man,ros$main)
      }

      consensus_rules<- nrls_man%>% distinct(features,levels,decision,.keep_all = T)
      #nrls<-mergeRBMs(rbms,defClass = eval(parse(text=paste("data$", decision_var, sep = ""))))
      return(consensus_rules)
    }else(return(ros))

  }else(return(ros))

}


##IMPORTANT####
#running consensus genetic model
#nrls<-genetic_consensus(data = data,num_of_iter = 2)
# saveRDS(nrls,'../../data/Rosetta_results/nrls_1.rds')








