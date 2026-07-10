#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

set_wd_to_script_dir <- function() {
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && !is.null(rstudioapi::getActiveDocumentContext()$path) &&
      rstudioapi::getActiveDocumentContext()$path != "") {
    path <- rstudioapi::getActiveDocumentContext()$path
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_flag <- "--file="
    match <- grep(file_flag, args)
    if (length(match) > 0) {
      path <- sub(file_flag, "", args[match])
    } else {
      stop("Could not determine script path. Run via Rscript or from an RStudio-opened file.")
    }
  }
  setwd(dirname(normalizePath(path)))
  message("Working directory set to: ", getwd())
}

set_wd_to_script_dir()
getwd()
set.seed(0)

library(org.Mmu.eg.db)


intersecting_colnames<-function(){
  #MCFS selection
  res<-readRDS('../../data/Mcfs_results/mcfs_D7_gp1.rds')
  res_1<-readRDS('../../data/Mcfs_results/mcfs_D133_gp1.rds')
  sig_genes<- res$RI$attribute[1:500]
  sig_genes_1<- res_1$RI$attribute[1:500]
  union(sig_genes,sig_genes_1)
  intersect(sig_genes,sig_genes_1)
  col<-intersect(sig_genes,sig_genes_1)
  return(col)
}


col<-intersecting_colnames()
col<-gsub('^X','ENSMMUG',col)
col<-gsub('[.]','-',col)

annot <- AnnotationDbi::select(org.Mmu.eg.db,
                               keys=col,
                               columns='ENSEMBL',
                               keytype="SYMBOL")

annot$ENSEMBL<- ifelse(is.na(annot$ENSEMBL), annot$SYMBOL,annot$ENSEMBL)
ensembl_ids<-annot$ENSEMBL

library(readxl)
DDE<-read_excel('../../data/Barrenäs et al/DDE_genes.xlsx')
DDE<-DDE$Ensembl.RM.ID

dde_intersect<- intersect(DDE, ensembl_ids)
print(paste0('Matching DDE gene is', dde_intersect,', SYMBOL:', annot$SYMBOL[annot$ENSEMBL==dde_intersect]))

DE<-read_excel('../../data/Barrenäs et al/DE_genes.xlsx')
DE<-DE$Ensembl.RM.ID

IL15<-read_excel('../../data/Barrenäs et al/IL15-response.xlsx')
IL15<-IL15$Ensembl.RM.ID



de_intersect<- intersect(DE, ensembl_ids)
for (i in de_intersect){
  #print(paste0('Matching DE gene is ', i,', SYMBOL:', annot$SYMBOL[annot$ENSEMBL==i]))
  print(annot$SYMBOL[annot$ENSEMBL==i])
}

print(annot$SYMBOL[annot$ENSEMBL %in% de_intersect])


il15_intersect<- intersect(IL15, ensembl_ids)
for (i in il15_intersect){
  #print(paste0('Matching DE gene is ', i,', SYMBOL:', annot$SYMBOL[annot$ENSEMBL==i]))
  print(annot$SYMBOL[annot$ENSEMBL==i])
}

intersect(annot$SYMBOL[annot$ENSEMBL %in% de_intersect],annot$SYMBOL[annot$ENSEMBL %in% il15_intersect])

print(annot$SYMBOL[annot$ENSEMBL %in% il15_intersect])





