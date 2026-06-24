#' @title boost_check_case
#' @author Girish Pulinkala
#' @description checks if the mcfs result is worth boosting
#' @param mcfs_result  A dataframe consiting of RI values from mcfs run.
#' @param col_size An integer value containng the nummber  in df
#' @return A integer specifying which case to run for boosting

boost_check_case<- function(mcfs_result,col_size){
  if(mcfs_result$cutoff_value < 10){
    message('Significant features from MCFS less than 10')
    message('\nBoosting Features will be done in order of significance from MCFS')
    if(col_size > 500){
      case_num<-2
      message('\nUpto 500 features')
    }else{case_num<-1
    message('\nUpto number of features present in data')}
  } else if(mcfs_result$cutoff_value > 10 & mcfs_result$cutoff_value < 100 ){
    message('Significant features from MCFS more than 10 but less than 100')
    message('\nBoosting Features will be done in order of significance from MCFS')
    if(col_size > 500){
      case_num<-4
    }else{case_num<-3
    message('\nUpto number of features present in data')}
  } else{
    message('Significant features from MCFS more than 100')
    message('\nBoosting Features will not be done')
    case_num<-5
  }
  return(case_num)
}

#' @title select_cols_boost
#' @author Girish Pulinkala
#' @description Selects the number of columns to boost. Since r.roseeta has a limit it selects a maximu of 500 columns
#' @param mcfs_result  A dataframe consiting of RI values from mcfs run.
#' @param col_size An integer value containng the nummber  in df
#' @param boost_case_num The output from the boost_check_case function
#' @return A vector of values containing columns that can be used to train your model.

select_cols_boost<-function(boost_case_num,col_size,mcfs_result){
  if(boost_case_num==1 | boost_case_num==3){
    cols<- mcfs_result$RI[1:(col_size-1),'attribute']
  }else if(boost_case_num==2 | boost_case_num==4){
    cols<- mcfs_result$RI[1:500,'attribute']
  }else(boost_case_num==5)( cols<- mcfs_result$RI[1:mcfs_result$cutoff_value,'attribute'])
  return(cols)
}

#' @title quality_check
#' @author Girish Pulinkala
#' @description does a quality check to find the accuracte or best model.
#' @param ros  A output from the rosetta function
#' @param quality A vector containing quality of all models from boosting
#' @return A boolean vaue inicating whether the current run was the best run.

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

#' @title quality_check
#' @author Girish Pulinkala
#' @description does a quality check to find the accuracte or best model.
#' @param ros  A output from the rosetta function
#' @param quality A vector containing quality of all models from boosting
#' @return A boolean vaue inicating whether the current run was the best run.


boost_ros<-function(ros_data,step=2,discrete=FALSE,
                    reducer='Johnson',
                    JohnsonParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95),
                    GeneticParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple")){
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
      ros<-rosetta(ros_data,reducerDiscernibility = 'Object',roc=TRUE,discrete = discrete,reducer = reducer,clroc=decision_roc,ruleFiltration = TRUE,JohnsonParam=JohnsonParam, GeneticParam=GeneticParam)
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
      ros<-rosetta(trim_ros_data,reducerDiscernibility = 'Object',roc=TRUE,discrete = discrete,reducer = reducer,clroc=decision_roc,ruleFiltration = TRUE,JohnsonParam=JohnsonParam, GeneticParam=GeneticParam)
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


#' @title mcfs_to_rosetta
#' @author Girish Pulinkala
#' @description Does feture boosting and builts rule-based models from mcfs input
#' @param ros  A output from the rosetta function
#' @param quality A vector containing quality of all models from boosting
#' @return A boolean vaue inicating whether the current run was the best run.


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


