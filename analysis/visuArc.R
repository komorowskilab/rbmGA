set.seed(0)
#Load libraries
#library(VisuNet)
if (!require(arcdiagram)) devtools::install_github('gastonstat/arcdiagram')
library(arcdiagram)


visuArc <- function(net, decision, feature=NULL,label=NULL){
  #load libraries if you only copy the function
  #library(VisuNet)
  # if (!require(arcdiagram)) install.packages('arcdiagram')
  # library(arcdiagram)
  
  #make df:  a data frame containing network connection in columns
 # df<- eval(parse(text=paste("net$", decision,'$nodes' ,sep = "")))
  df<-net[[decision]][['nodes']]
  cols_select<- c('label','DiscState','NodeConnection')
  df<- subset(df,select = cols_select)
  rownames(df)  <-make.names(df[,'label'], unique = TRUE)
  df <- df[,-1]
  
  #scale Node connection fro 0 to 1
  scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
  df$NodeConnection <- scale_values(df$NodeConnection)
  
  
  #check for feature
  if(is.null(feature)){
    feature <- rownames(df[which.max(df[,'NodeConnection']),])
  }
  
  discLev<-df[feature,'DiscState']
  feature2 <- gsub("\\s*\\([^\\)]+\\)","",as.character(feature))
  edges <- net[[decision]]$edges[,1:4]
  edges_from <- as.character(edges$from)
  edges_to <- as.character(edges$to) ## gsub("\\=.*","",as.character(net$short$edges$to))
  indx <- unique(c(which(edges_from == paste0(feature2,"=",discLev)),which(edges_to == paste0(feature2,"=",discLev))))
  edges2 <- edges[indx,]
  M <- as.matrix(edges2[,1:2])
  
  #colors
  valsCols <- round(edges2$connNorm, digits = 1)*10
  edgesCon <- as.numeric(as.factor(as.character(valsCols)))
  colors <- colorRampPalette(c("gainsboro","lavender","darkorchid3"))(length(unique(valsCols)))[edgesCon]
  labelsNodes <- unique(c(t(M)))
  colNodes <- c("#56B4E9","#999999","#E69F00")[as.numeric(sub(".*=","",labelsNodes))] #EDIT THIS WITH YOUR COLOR CODING
  
  #nodes values
  nodesDlev<- df[,'DiscState']
  nodesNams <- gsub("\\s*\\([^\\)]+\\)","",as.character(rownames(df)))
  nodes2 <- paste0(nodesNams,"=",nodesDlev)
  nodeSize <- round(df[match(labelsNodes, nodes2),2],digits = 1)*10
  
  ordV <- c(1,order(edgesCon, decreasing = T)+1)
  
  colsLabs <- rep("black",length(labelsNodes))
  colsLabs[which(labelsNodes == paste0(feature2,"=",discLev))] <- "red"
  
  if(is.null(label)){
    label=decision
  }
  
  arcplot(M, lwd.arcs=edgesCon, col.arcs = colors, col.nodes = colNodes, labels=sub("=.*","",labelsNodes),
          ordering=ordV, col.labels=colsLabs, cex.labels=0.7, font=1, lwd.nodes = nodeSize, horizontal = FALSE,main = label,adj=1)
}

#Manual
#decision : decision
#net : Visunet output
#feature : central or interested hub feature. P.S. It should have connections!!!!!!! If null it chooses the central hub by default.

##Program 
#ros<-rosetta(autcon)
#vis<-visunet(ros$main)
#visuArc(vis,'control',feature='PPOX')
