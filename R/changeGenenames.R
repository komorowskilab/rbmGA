#' @title chanege gene names
#' @author Girish Pulinkala
#' @description converts ensembl IDs to symbols by preserving the ensembl ids that do not have any Symbol yet.
#' @param data decision table
#' @param transpose A boolean value that indicates whether to transpose the data table
#' @param ensembl_id_col A integer value denotng the column where the column is present.
#' @param remove_cols A list/ vector that contains index values of columns to be removed in the process.
#' @examples
#' #change_gene_names(data,transpose=TRUE,ensembl_id_col=1,remove_cols=c(2))
#'


change_gene_names<- function(data,transpose=TRUE,ensembl_id_col=1,remove_cols=c(2)){


  if(transpose==TRUE){
    data<-as.data.frame(t(data))
    Gene.stable.ID <- rownames(data[ensembl_id_col,])
  }

  Gene.stable.ID<-names(data)[ensembl_id_col]

  #data<- data %>% rownames_to_column(.,var='ENSEMBL')  #Optional
  library(org.Mmu.eg.db)
  mapped_df <- data.frame(
    "entrez_id" = mapIds(
      # Replace with annotation package for the organism relevant to your data
      org.Mmu.eg.db,
      keys = data[[ensembl_id_col]],
      # Replace with the type of gene identifiers in your data
      keytype = "ENSEMBL",
      # Replace with the type of gene identifiers you would like to map to
      column = "SYMBOL",
      # This will keep only the first mapped value for each Ensembl ID
      multiVals = "first"
    )
  ) %>%
    # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
    # drop that from the data frame
    #dplyr::filter(!is.na(entrez_id)) %>%
    # Make an `Ensembl` column to store the row names
    tibble::rownames_to_column("Ensembl") %>%
    # Now let's join the rest of the expression data
    dplyr::inner_join(data, by = c("Ensembl" = Gene.stable.ID))

  mapped_df$entrez_id <- ifelse(is.na(mapped_df$entrez_id), substring(mapped_df$Ensembl,8,nchar(mapped_df$Ensembl)),mapped_df$entrez_id)
  #mapped_df<-subset(mapped_df, select= -c(Ensembl))
  mapped_df<-subset(mapped_df, select= -c(Ensembl,HGNC.symbol))
  rownames(mapped_df) = make.names(mapped_df$entrez_id, unique=TRUE)
  mapped_df<-subset(mapped_df, select= -c(entrez_id))
  mapped_df<-as.data.frame(t(mapped_df))
  return(mapped_df)
}

change_gene_names(count_data,ensembl_id_col = 1,transpose = F)
