#' @title Step Wise Selection Method
#' @author Girish Pulinkala
#' @description Does step-wise selection to select the best and accurate model
#' @import R.ROSETTA
#' @import data.table
#' @import tibble
#' @import stringr
#' @import grid
#' @import dplyr
#' @import ggplot2
#'
#' @param data decision table used in the training
#' @param iterations A rule set that was generated using the training set
#' @param reducer A character containing name of reducer method: Johnson or Genetic. Default is Johnson.
#' @param discrete Logical. Set TRUE for discrete data. Default is FALSE.
#' @param discreteMethod A character containing discretization method: EqualFrequency, MDL, Naive, SemiNaive or BROrthogonal. Default is EqualFrequency.
#' @param discreteParam A vector containing discretization parameters. May be of different length and values. See examples.
#' @param discreteMask Logical. Set FALSE to disable discretization mask. Default is TRUE.
#' @param cvNum A numeric value of the cross-validation number. Default is 10.
#' @param JohnsonParam A vector containing Johnson reducer parameters.
#' @param GeneticParam A vector containing Genetic reducer parameters. description
#'
#' @examples
#' # example code
#' # df<-fread('../path/to/decision_table)
#' # result<- step_wise(df, reducer='Johnson,discrete=FALSE)
#' # print(result[[2]])


step_wise<-function(data,step=5,reducer = 'Johnson',discrete=FALSE,discreteMethod = 'EqualFrequency',
                    cvNum=10,ruleFiltration = TRUE,discreteParam =3,
                    JohnsonParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95),
                    GeneticParam = list(Modulo=TRUE, BRT=FALSE, BRTprec=0.9, Precompute=FALSE, Approximate=TRUE, Fraction=0.95, Algorithm="Simple")){

  #foundational function
  quality_check<- function(ros,quality){
    if(length(quality)==1){
      return (TRUE)
    }else if(ros$quality$accuracyMean > max(quality)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }


  #Main script
  col_size<-ncol(data)-1
  decision_var<-names(data)[ncol(data)]
  decision_roc<-names(table(eval(parse(text=paste("data$", decision_var, sep = ""))))[1])

  df<-data.frame()

  #print(col_size)
  count=10
  while(count<=plyr::round_any(col_size, step, f = ceiling) ){
    if(count>=col_size){
      ros<-rosetta(data,reducerDiscernibility = 'Object',roc=TRUE,clroc=decision_roc,fallBackClass = decision_roc,discrete = discrete,discreteMethod = discreteMethod,discreteParam = discreteParam,reducer = reducer,ruleFiltration = ruleFiltration,cvNum = cvNum)
      quality<-append(quality,ros$quality$accuracyMean)

      ######################################
      if(quality_check(ros,quality)==TRUE){
        ros_max<-ros
        data_max <-data;
        recal<-recalculateRules(data,ros$main,discrete = discrete)
      }
      ######################################

      response_var <- names(data)[ncol(data)]
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
      trim_data<-data[,1:count]
      #print(dim(trim_data))
      trim_data$decision<-eval(parse(text=paste("data$", decision_var, sep = "")))
      trim_data$decision<-factor(trim_data$decision)
      ros<-rosetta(trim_data,reducerDiscernibility = 'Object',roc=TRUE,clroc=decision_roc,fallBackClass = decision_roc,discrete = discrete,discreteMethod = discreteMethod,discreteParam = discreteParam,reducer = reducer,ruleFiltration = ruleFiltration,cvNum = cvNum)
      quality<-append(quality,ros$quality$accuracyMean)


      ######################################
      if(quality_check(ros,quality)==TRUE){
        ros_max<-ros
        data_max <-trim_data
        recal<-recalculateRules(trim_data,ros$main,discrete = discrete)
      }
      #####################################

      #metrics
      df[count,'Nof']<-count
      df[count,'RBM_accuracy']<-ros$quality$accuracyMean
      df[count,'RBM_accuracy_std']<-ros$quality$ROC.AUC.MEAN

      #Counter
      count<-count+step

    }

  }
  step_wise_plot<-na.omit(as.data.frame(df)) %>% ggplot()  +
    geom_point(aes(Nof,RBM_accuracy,shape='meanAccuracy',color='meanAccuracy'),size=5,,show.legend = T)+
    geom_point(aes(Nof,RBM_accuracy_std,shape='meanROC',color='meanROC'),size=5,show.legend = T) +
    scale_x_continuous(breaks = (as.data.frame(df)$Nof))+   scale_y_continuous(breaks = seq(0,1,0.1)) + expand_limits(x = 10, y = 0) +
    scale_shape_manual(values = c("meanAccuracy" = 3, "meanROC" = 4)) +
    scale_color_manual(values = c("meanAccuracy" = "black", "meanROC" = "red")) +
    labs(x='Number of top features selected ', y='Value', shape='Metric', color='Metric')+
    ggtitle('Step-wise Selection_Model') + #EDIT THIS To give title
    theme(axis.text.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.text.y = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.x = element_text(family = 'sans', face='plain', colour='black', size=12),axis.title.y = element_text(family = 'sans', face='plain', colour='black', size=12),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme_minimal()

  print(step_wise_plot)

  return(list(data_max,ros_max,recal,step_wise_plot))
}



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



