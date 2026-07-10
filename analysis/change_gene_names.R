#REQUIRES NORM DATA AS GIVEN BY GALE LAB
data_prep<- function(data,transpose=TRUE){


  if(transpose==TRUE){
    data<-as.data.frame(t(data))

  }


  #data<- data %>% rownames_to_column(.,var='ENSEMBL')  #Optional
  library(org.Mmu.eg.db)
  mapped_df <- data.frame(
    "entrez_id" = mapIds(
      # Replace with annotation package for the organism relevant to your data
      org.Mmu.eg.db,
      keys = data$Gene.stable.ID,
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
    dplyr::inner_join(data, by = c("Ensembl" = "Gene.stable.ID"))

  mapped_df$entrez_id <- ifelse(is.na(mapped_df$entrez_id), substring(mapped_df$Ensembl,8,nchar(mapped_df$Ensembl)),mapped_df$entrez_id)
  #mapped_df<-subset(mapped_df, select= -c(Ensembl))
  mapped_df<-subset(mapped_df, select= -c(Ensembl,HGNC.symbol))
  rownames(mapped_df) = make.names(mapped_df$entrez_id, unique=TRUE)
  mapped_df<-subset(mapped_df, select= -c(entrez_id))
  mapped_df<-as.data.frame(t(mapped_df))
  return(mapped_df)

}

