#' @title create_indexed_file
#' @author Girish Pulinkala
#' @description creates indexed log file if a log file exists
#' @param filename A input from check_log_file
#' @param extension A input from check_log_file
#' @param file_path A input from check_log_file

create_indexed_file = function(filename,extension='.txt', file_path){

  if(!file.exists(paste0(file_path,filename,extension))){return(paste0(filename,extension))}
  i=1
  repeat {
    f = paste0(filename,i,extension)
    if(!file.exists(paste0(file_path,f))){return(f)}
    i=i+1
  }
}



#' @title check_log_file
#' @author Girish Pulinkala
#' @description Checks if a log file exists in the folder
#' @param mcfs_out_path  A file path to save mcfs files

check_log_file<-function(mcfs_out_path){

  if(file.exists(mcfs_out_path)){
    cat('Creating new out file since out.txt already exists')
    new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'/Mcfs_results/'))
    file.create(paste0(getwd(),'/Mcfs_results/',new_indexed_file))
    return(paste0(mcfs_out_path,'/',new_indexed_file))
  }
  else{
    new_indexed_file<-create_indexed_file('out',extension = '.txt',file_path= paste0(getwd(),'/Mcfs_results/'))
    file.create(paste0(getwd(),'/Mcfs_results/',new_indexed_file))
    return(paste0(mcfs_out_path,'/',new_indexed_file))
  }

}


#' @title runMCFS
#' @author Girish Pulinkala
#' @description Runs Monte Carlo Feature selection and saves outputs and logs
#' @param mcfs_df A dataframe or a decision table with decision as the last column.
#' @param out_name A string filename for the final mcfs.rds file, for e.g. 'autcon'
#' @param permutations A integer value which specfies the number of permutations It is from the original mcfs function.
#' @param set_response_var A boolean value, which says if the user wamts to specifiy the decision column. by default it is set to FALSE. Use this if the your last column is not the decision column
#' @param response_var A string value to specify the decision column name. By default is set to FALSE and will automatically set the last column to be decision column
#' @param plot_res A boolean value that sepcifies whether to plot the graphs from the mcfs runs. By deafult it is set to TRUE.
#' @param projections An integer value specifying number of projections. Refer to rmcfs documentation.
#' @param projectionSize An integer value specifying number of projectionSize. Refer to rmcfs documentation.
#' @param balance An integer value to specify undersampling. Refer to rmcfs documentation.
#' @export
#' @import data.table
#' @import dplyr
#' @import rJava
#' @import rmcfs
#' @importFrom grDevices dev.off
#' @importFrom stats as.formula
#' @importFrom utils capture.output
#' @examples
#' # example code
#' # dt <- read.delim('alizadeh.csv)
#' # runMCFS(dt, 'mcfs_run1',permutations=20,set_response_var=FALSE,response_var=NULL,plot_res=TRUE)

runMCFS<-function(mcfs_df,out_name,permutations=3,set_response_var=FALSE,response_var=NULL,plot_res=TRUE,projections='auto', projectionSize='auto',balance=2){

  # requires create_indexed_file
  # requires check_log_file
  ## mcfs_out_path : file path for saving mcfs file
  ## plot_out_path : file path for saving mcfs plots and graphs

  #Directories and files
  dir.create(file.path(getwd(), 'Mcfs_results'),showWarnings = FALSE)
  mcfs_out_path<- file.path(getwd(), 'Mcfs_Results')


  mcfs_df<-as.data.frame(mcfs_df)
  cat("Checking input arguments...\n")
  if(is_empty(out_name)){
    stop('output name is empty')
  }


  cat('Checking for out file')
  out_file<-check_log_file(mcfs_out_path)

  #fixdata for mcfs
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
    formula <- stats::as.formula(paste(response_var, "~ ."))
  }

  cat(paste0('The Response varibale is : ',response_var))
  cat(paste0('\nStarted running Monte Carlo Feature Selection at ',Sys.time()))

  #run monte carlo feature
  mcfs_result <- rmcfs::mcfs(formula, mcfs_df, projections=projections, projectionSize=projectionSize, splits=5, cutoffPermutations=permutations,threadsNumber= 8,balance=balance)

  #Save mcfs in R data format
  saveRDS(mcfs_result,paste0(mcfs_out_path,'/',out_name,'.rds'))

  #Plot Graphs and results
  if(plot_res==TRUE){
    dir.create(file.path(getwd(), 'Mcfs_results/plots'),showWarnings = FALSE)
    plot_out_path<-file.path(getwd(), 'Mcfs_results/plots')
    gid <- rmcfs::build.idgraph(mcfs_result, size = 6, size_ID = 12, orphan_nodes = TRUE)
    rmcfs::export.plots(mcfs_result, mcfs_df, idgraph = gid, path = plot_out_path, label = "mcfs_plot_", color = "darkgreen",plot_format = 'png')
    grDevices::dev.off() #Not sure if this is the right place but it works
  }

  cat(paste0('Mcfs file Saved in Directory : ',paste0(mcfs_out_path,'/',out_name,'.rds')))
  cat(paste0('\nMCFS run exited with response variable: ',response_var))
  cat('\nCheck Log file in directory for execution IFNO')
  cat('######## Monte Carlo Feature Selection ############',file=out_file,sep="\n",append=TRUE)
  cat(paste0('Files and plots saved in Directory : ',mcfs_out_path),file=out_file,sep="\n",append=TRUE)
  cat(utils::capture.output(print.mcfs(mcfs_result)),file = out_file,sep="\n",append=TRUE)
  #cat(capture.output(),file = out_file)
}

