#' @title TMMnormalization
#' @author Girish Pulinkala
#' @description TMM normalisation for decision system
#' @param dt  A decision table
#' @param Groups A string value denoting the condition or treatment or classification
#' @returns A dataframe for ML
#' @export
#' @import edgeR
#' @import sva
#' @examples
#' # norm_dt<- TMMnormalisation(dt)


TMMnormalisation <- function(dt, Groups='attr_oncotree_disease_code'){

  counts<- dt %>% dplyr::select(-Groups) %>% t()
  group <- dt %>% dplyr::select(Groups)

  d <- DGEList(counts=counts,group = group$Groups)

  TMM<- calcNormFactors(d,method='TMM')
  # print(plotMDS(TMM, method="bcv", col=as.numeric(d$samples$group)))
  # print(legend("bottomleft", c('iALL','pALL'), col=1:3, pch=20))

  TMM <- cpm(TMM, log = FALSE, normalized.lib.sizes=TRUE)
  mod<-model.matrix(~group[[Groups]],data=d)
  TMM <- TMM[rowMeans(TMM) != 0,]
  svaresult<-sva(TMM,mod)
  bat<-removeBatchEffect(TMM,sva=svaresult)
  bat<-as.data.frame(t(bat))
  rownames(bat)<-rownames(group)

  if(all(rownames(bat) == rownames(groups))){
    bat$Group <- group[[Groups]]
  }else{
    message('The rownames were not in order')
    bat <- bat[rownames(group),]
    bat$Group <- group[[Groups]]
  }
  return(bat)

}
