#' @title Rule Characterization
#' @author Girish Pulinkala
#' @description Characterize the rules into the hubs.
#' @import pheatmap
#' @import data.table
#' @import tibble
#' @import stringr
#' @import dplyr
#' @import R.ROSETTA
#'
#' @export
#'
#' @param df decision table used in the training
#' @param mod_recal recalutated rules from R.ROSETTA
#' @param group
#' @param type To descibe the hub clusters identified in the analysis
#' @param TP timepoint
#' @examples
#' # example code
#' #rule_characterization(ros_data,mod_recal,group=groups[i],type='Rules',TP=TP)
#'

rule_characterization<-function(df,recal,group=NULL,type='Cluster',TP='D7'){
  library(pheatmap)
  #process data for clustering, i.e. making a matrix with binary values 0 and 1
  dataM <- data.frame(matrix(ncol = length(recal$features[1:nrow(recal)]), nrow = length(row.names(df))))
  rownames(dataM)<-rownames(df)
  colnames(dataM)<-recal$features
  #dataM<-rbind(recal$C,dataM)

  #assigning 1 to objects satisfying rules
  for(i in 1:length(recal$features[1:nrow(recal)])){
    list_of_features<-unlist(strsplit(recal$supportSetLHS[i],split=','))
    if(!is_empty(list_of_features)){
      for(j in 1:length(list_of_features)){
        if(list_of_features[j] %in% rownames(dataM))(dataM[list_of_features[j],i]<-recal$accuracyRHS[i])
        # if(list_of_features[j] %in% rownames(dataM))(dataM[list_of_features[j],i]<-1)
      }
    }
  }

  #replacing not satisfying rules as 0
  dataM[is.na(dataM)] <- 0

  col_anno<-setNames(recal[[type]], recal$features)
  summary <- sapply(unique(col_anno), function(cat) {
    cols <- names(col_anno[col_anno == cat])
    row_totals <- rowSums(dataM[, cols, drop=FALSE])

    cat_total <- row_totals/length(cols)

    return((cat_total))

  })

  print(summary)

  return(plotting(summary,group=group,type=type,TP=TP))

}

plotting<-function(dataM,group=NULL,type='Cluster',TP=TP){
  gender='Male'
  group=TP


  cluster_cols<-c('Cell_Death'='royalblue','Inf_Control'='red','Other_in_Prot'='skyblue','NonProt'='yellow2')

  reshape2::melt(dataM) %>% filter(if_all(everything(), ~. != 0)) %>% ggplot(aes(as.character(Var1), value, fill = Var2)) +
    geom_bar(stat = "identity") +
    coord_flip() + theme_bw()


  p1 <- dataM %>%
    reshape2::melt() %>%
    left_join(
      ros_data %>%
        mutate(Var1 = as.numeric(rownames(.))) %>%
        dplyr::select(Var1, protectionStatus),
      by = 'Var1'
    ) %>%
    filter(protectionStatus == 'Prot') %>%
    mutate(Var2 = factor(Var2, levels = rev(names(cluster_cols)))) %>%
    ggplot(aes(x = as.character(Var1), y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.1) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
    scale_fill_manual(values = cluster_cols) +
    coord_flip() +
    labs(
      y = NULL,
      x = NULL,
      fill = type
    ) +
    ggtitle(paste0(gender, ' Protected Animals_', group)) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )


  p2 <- dataM %>%
    reshape2::melt() %>%
    left_join(
      ros_data %>%
        mutate(Var1 = as.numeric(rownames(.))) %>%
        dplyr::select(Var1, protectionStatus),
      by = 'Var1'
    ) %>%
    filter(protectionStatus == 'NonProt') %>%
    mutate(Var2 = factor(Var2, levels = rev(names(cluster_cols)))) %>%
    ggplot(aes(x = as.character(Var1), y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.1) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
    scale_fill_manual(values = cluster_cols) +
    coord_flip() +
    labs(
      y = NULL,
      x = NULL,
      fill = type
    ) +
    ggtitle(paste0(gender, ' Non-Protected Animals_', group)) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )


  return(list(p1,p2))
}







