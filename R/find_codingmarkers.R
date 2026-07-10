#' @title findCodingmarkers
#' @author Girish Pulinkala
#' @description Find protein coding and trim the dataframe from all other or noncoding regions
#' @param dt  A data frame containing decision table. The last column of decision should be missing at this point. For experts, this is an information system.
#' @param keytype the keytype that matches the keys used. For the select methods, this is used to indicate the kind of ID being used with the keys argument. For the keys method this is used to indicate which kind of keys are desired from keys
#' @returns A dataframe with only protein coding regions, if transposed can be used for differnetial expressions too
#' @export
#' @import dplyr
#' @import tibble
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @examples
#' # path <- c("../data/data_all_old/")
#' # dt<- readFeaturecounts(path,makedt=TRUE)
#' # dt<-findCodingmarkers(dt,keytype="SYMBOL")


findCodingmarkers<-function(dt,keytype="SYMBOL"){
  keytypeCols <- c("SYMBOL",'GENETYPE')

  if(keytype!="SYMBOL"){
    keytypeCols <-append(keytype, keytypeCols)

  }

  annot <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys=colnames(dt),
                                 columns=keytypeCols,
                                 keytype=keytype)

  gene<- annot[which(annot$GENETYPE == 'protein-coding'),]

  dt<-subset(dt,select=gene$SYMBOL)

  return(dt)
}



#' @title removeRPgenes
#' @author Girish Pulinkala
#' @description Remove ribosomal protein coding genes from the decision table based on a regular expression value.
#' @param dt  A data frame containing decision table. The last column of decision should be missing at this point. For experts, this is an information system.
#' @returns A dataframe with no Ribosomal protein coding genes.
#' @export
#' @import stringr
#' @examples
#' # path <- c("../data/data_all_old/")
#' # dt<- readFeaturecounts(path,makedt=TRUE)
#' # dt<-findCodingmarkers(dt,keytype="SYMBOL")
#' # dt <- removeRPgenes(dt)


removeRPgenes <- function(dt){

  RPL_pattern<-'^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA'

  RPL_genes<-subset(dt,select=(stringr::str_detect(names(dt),RPL_pattern))) %>% names()
  message('Removed Ribosomal coding genes:')
  message(paste(RPL_genes,collapse=','))
  dt<-dt %>% dplyr::select((-RPL_genes))

  return(dt)

}


#' @title addDecision
#' @author Girish Pulinkala
#' @description Remove non-coding, pseudo, rRNAs, scRNA, snoRNA, snRNA, unknown and other genes based on org.Hs.eg.db
#' @param dt  A data frame containing decision table. The last column of decision should be missing at this point. For experts, this is an information system.
#' @param metadata_path A string pointing to the metadta file in your system. for e.g. '../data/SAMPLE_INFO.txt'
#' @param trim A boolean value to be set to TRUE by people using St Judes dataset. By default is set to FALSE
#' @param sample_column A integer value containing the column of the sample names in the metadata
#' @param target_column A integer value containing the column of the decision value to be appended to the decision table.
#' @returns A dataframe which is a true decision table with decision in the last column and could be used in ML applications.
#' @export
#' @import data.table
#' @import dplyr
#' @import org.Hs.eg.db
#' @examples
#' # dt<- addDecision(dt, metadata_path="../data/SAMPLE_INFO.txt",
#' # trim=FALSE, sample_column=4 ,target_column=23)




addDecision<-function(dt, metadata_path="../data/SAMPLE_INFO.txt", trim=FALSE, sample_column=4 ,target_column=23){

  metadata<-fread(metadata_path)

  sample_name<-names(metadata)[sample_column]

  #TRIMMING IS APPLICABLE TO ONLY ST JUDES DATASET.
  #MAKE SURE YOUR COLUM TO JOIN IS SAME AS YOUR ROWNAMES
  if(trim==TRUE){
    trim_meta <- function (x) gsub('.RNA-Seq.feature-counts','', x)
    dt[['sample_name']]<- trim_meta(rownames(dt))
  }

  dt <-dt %>% dplyr::left_join(metadata %>% dplyr::select(sample_column,target_column),by=c('sample_name'=sample_name)) %>% dplyr::distinct() %>% tibble::column_to_rownames(.,'sample_name')

  return(dt)

}







