# rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #will throw an error when not file not saved
# getwd()
# set.seed(0)

# if (!require(package)) install.packages('package')
# library(package)


library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(plyr)
#library(org.Hs.eg.db)
library(R.ROSETTA)
library(MVN)
library(ggplot2)
library(stringr)
library(pheatmap)
library(funModeling)

library(gridExtra)
library(grid)
library(ComplexHeatmap)


#start_time<-Sys.time()
cluster_rules<-function(training_df,recal,show_colnames=FALSE,show_rownames=FALSE,support=5,flag=F){
  library(pheatmap)
  #process data for clustering, i.e. making a matrix with binary values 0 and 1
  dataM <- data.frame(matrix(ncol = length(recal$features[1:nrow(recal)]), nrow = length(row.names(training_df))))
  rownames(dataM)<-rownames(training_df)
  #colnames(dataM)<-paste(recal$features,recal$levels,sep='_')
  for(i in 1:nrow(recal)){
    names(dataM)[i]<- paste(unlist(paste(unlist(strsplit(as.character(recal[i,]$features),",")) ,unlist(strsplit(as.character(recal[i,]$levels),",")),sep='=')),collapse=',')
  }
  
  #assigning 1 to objects satisfying rules
  for(i in 1:length(recal$features)){
    list_of_objects<-unlist(strsplit(recal$supportSetLHS[i],split=','))
    for(j in 1:length(list_of_objects)){
      if(list_of_objects[j] %in% rownames(dataM))(dataM[list_of_objects[j],i]<-1)
    }
  }
  
  #replacing not satisfying rules as 0
  dataM[is.na(dataM)] <- 0
  
  #decision variable of the dataset
  decision_var<-names(training_df)[ncol(training_df)]
  ann <- data.frame( eval(parse(text=paste("training_df$", decision_var, sep = ""))))
  colnames(ann) <- 'Decision'
  rownames(ann)<- rownames(dataM)
  newdf<-dataM[rowSums(dataM[])>1,colSums(dataM[])>support]
  a <- filter(ann, rownames(ann) %in% rownames(newdf))
  annoCol <- list(category = unique(ann$Decision))
  
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  p<- pheatmap::pheatmap(as.matrix(newdf), annotation_row=a, fontsize_col = 5 ,fontsize_row = 10, border_color = 'white',annotation_colors = annoCol,cluster_cols = TRUE,show_rownames = show_rownames, show_colnames = show_colnames, cluster_rows = TRUE,color = c('grey88','gray39'),cutree_rows = 4,cutree_cols = 2,legend_breaks = c(0,1))
  print(pheatmap::pheatmap(as.matrix(newdf), annotation_row=a, fontsize_col = 5 ,fontsize_row = 10, border_color = 'white',annotation_colors = annoCol,cluster_cols = TRUE,show_rownames = show_rownames, show_colnames = show_colnames, cluster_rows = TRUE,color = c('grey88','gray39'),cutree_rows = 4,cutree_cols = 2,legend_breaks = c(0,1)))
  setHook("grid.newpage", NULL, "replace")
  grid.text("Rules", y=-0.07, gp=gpar(fontsize=15))
  grid.text("Samples", x=-0.07, rot=90, gp=gpar(fontsize=15))
  
  
  
  
  
  annoCol <- list(Decision = c("NonProt" = "#999999", "Prot" = "#E69F00"))
  print(Heatmap(as.matrix(newdf) ,right_annotation = rowAnnotation(Decision=a$Decision,col=annoCol),col=colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100),name='Value',column_names_rot = 45, show_column_names = F,column_km  = 2, row_km  = 2,heatmap_legend_param = list(title = "value", color_bar = "discrete",legend_gp = gpar(fill = 0:1))))
  
  
  if(flag==TRUE){
    annorow <- list(Cluster = c("Cell_death" = "lightblue", "Inf_Control" = "darkgreen"))
   # recal$C)
    #p<-Heatmap(as.matrix(newdf),top_annotation = HeatmapAnnotation(Cluster=recal$C,col=annorow),right_annotation = rowAnnotation(Decision=a$Decision,col=annoCol),col=colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100),name='Value',column_names_rot = 45, show_column_names = F,column_km  = 2, row_km  = 2,heatmap_legend_param = list(title = "value", color_bar = "discrete",legend_gp = gpar(fill = 0:1)))
    print(Heatmap(as.matrix(newdf),top_annotation = HeatmapAnnotation(Cluster=recal$C,col=annorow),right_annotation = rowAnnotation(Decision=a$Decision,col=annoCol),col=colorRampPalette((brewer.pal(n = 7, name="OrRd")))(100),name='Value',column_names_rot = 45, show_column_names = F,column_km  = 2, row_km  = 2,heatmap_legend_param = list(title = "value", color_bar = "discrete",legend_gp = gpar(fill = 0:1))))
  }
  
  
  
  #for supplementary data, 
  #return(p)
  
}

#cluster_rules(ros_data,recal)
#EXAMPLE
# ros<-rosetta(autcon)
# recal<-recalculateRules(autcon,ros$main)
# cluster_rules(autcon,recal,show_colnames = FALSE)





