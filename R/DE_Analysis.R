#' @title DEAnalysis from a De cision table
#' @author Girish Pulinkala
#' @description generates DE analysis for GE data
#' @param dt A Decision table
#' @param keep_val A integer  which specifies the rowMeans value
#' @param compare A string from the binary outcome that you want to compare with
#' @export
#' @import DESeq2
#' @import ggrepel


DEAnalysis <- function (dt,keep_val=20, compare='iALL'){

  names(dt)[ncol(dt)] <- 'Groups'

  counts<- dt %>% dplyr::select(-Groups) %>% t()
  group <- dt %>% dplyr::select(Groups)

  group$Groups <- sort(group$Groups)

  if(!all(colnames(counts) == rownames(groups))){
    counts <- counts[,rownames(group)]
  }

  names(group)
  counts_filtered <- counts[, rownames(group)]

  #DESEQ" run
  dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                colData = group,
                                design = ~Groups)

  keep <- rowMeans(counts(dds)) >= keep_val
  dds <- dds[keep,]

  other<-setdiff(group$Groups,compare)

  dds$Groups <- factor(dds$Groups, levels = c(other,compare))

  dds <- DESeq(dds)

  res <- results(dds)
  message(resultsNames(dds))

  coef <- resultsNames(dds)[2]

  #shrinking the LFCs
  resLFC <- lfcShrink(dds, coef=coef, type="apeglm")

  #Plot MA plot
  print(plotMA(res))
  print(plotMA(resLFC))

  print(plotCounts(dds, gene=which.min(res$padj), intgroup="Groups"))

  #significant genes
  resOrdered <- res[order(res$pvalue),]
  resSig <- subset(resOrdered, padj < 0.01)


  # Volcano plot
  df <- as.data.frame(res)

  df$Significance <- "NoSig"
  df$Significance[df$log2FoldChange > 1 & df$padj < 0.05] <- "Up"
  df$Significance[df$log2FoldChange < -1 & df$padj < 0.05] <- "Down"


  print(df %>% rownames_to_column(.,'gene_symbol')  %>% ggplot(aes(x=log2FoldChange, y=-log10(padj), color=Significance)) +
    geom_point(alpha=0.4) +
    geom_text_repel(data = filter(df %>% rownames_to_column(.,'gene_symbol') , padj < 0.05 & (abs(log2FoldChange) > 3)),
                    aes(label = gene_symbol),
                    size = 3,
                    show.legend = FALSE) +
    theme_minimal() +
    scale_color_manual(values=c("Up" = "#00AFBB",
                                "Down" = "red",
                                "No" = "grey"))+
    scale_x_continuous(breaks = seq(-10, 10, 2))) #+   ggtitle('DE genes in', compare,'vs',other))

    return(list(res,resSig))

}
