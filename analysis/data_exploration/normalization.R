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
########## -- Libraries --##########

library(limma)
library(edgeR)
library(stats)
library(factoextra)
library(umap)
library(Rtsne)
library(ggplot2)
library(ExpressionNormalizationWorkflow)
library(stringr)
library(pvca)
library(mclust)
library(fossil)
library(amap)
library(rrcov)
library(svglite)
library(data.table)
library(circlize)

# TITLE: Watkins_07
# AUTHOR: LEANNE WHITMORE ed. KS

########## -- Files & Objects --##########

# rhesus2human <- read.csv(
#   file = "./RNAseq analysis files/rhesus2human109v2_ks.csv",
#   header = TRUE,
#   stringsAsFactors = FALSE
# )

########## -- FUNCTIONS --##########

theme_Publication <- function(base_size = 14, base_family = "arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2), hjust = 0.5
      ),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = unit(0.6, "cm"),
      legend.margin = unit(0, "cm"),
      legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    ))
}

generate_folder <- function(foldername) {
  workDir <- getwd()
  subDir <- foldername
  results_path <- file.path(workDir, subDir)
  if (file.exists(subDir)) {
  } else {
    dir.create(results_path)
  }
  return(results_path)
}

vizualize_DE_genes_bp <- function(results, plot_file) {
  print("STATUS: Generating bar plot of number of DE genes...")
  results_t <- t(summary(results))
  results_t <- results_t[, -2]

  for (i in 1:(length(row.names(results_t)))) {
    results_t[i, 1] <- results_t[i, 1] * -1
  }

  DE <- as.data.frame(results_t)
  DE <- setnames(DE,
                 old = c("Var1", "Var2", "Freq"),
                 new = c("Timepoint", "grodwn", "DE_genes")
  )

  # Create plot
  ggplot(DE, aes(
    x = Timepoint, y = DE_genes, fill = grodwn,
    label = DE$DE_genes
  )) +
    geom_bar(stat = "identity", position = "identity") +
    # geom_text(size = 5, position = position_stack(vjust = 0) )+
    # theme_light() +
    theme_minimal() +
    scale_fill_manual(values = c("#0808c4", "#da9618")) +
    # xlab("Time point")
    ylab("Number of Differentially Expressed Genes") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 15)
    )
  ggsave(plot_file, width = 6, height = 4, units="in", dpi = 300, bg="white")
}

generate_boxplots_voom <- function(data, labels, filename, figres, maintitle, ylabtitle) {
  png(filename, width = 10, height = 8, units = "in", res = figres)
  # par(mar=c(1,1,1,1))
  minvalue <- min(data)
  maxvalue <- max(data)
  boxplot(data,
          labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
          ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
          frame = FALSE
  )
  dev.off()
}

generate_density_plot <- function(data, labels, filename, figres) {
  png(filename, res = figres)
  par(xpd = TRUE)
  if (length(labels) > 10) {
    plotDensities(data, legend = FALSE)
  } else {
    plotDensities(data,
                  legend = "topright",
                  inset = c(-0.2, 0), levels(labels)
    )
  }
  dev.off()
}

normalize_data <- function(CM2, targetfile, grodwn, control = FALSE) {

  # order target and count matrix so they are the same (THIS IS IMPORTANT)
  CM2 <- CM2[, rownames(targetfile)]

  # CHECK IF ORDER IS THE SAME
  if (all.equal(colnames(CM2), rownames(targetfile)) != TRUE) {
    print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
    print(rownames(targetfile))
    print(colnames(CM2))
  }
  # biological reps

  ##CHANGE FOR DIFFERENT ANALYSES
  #MamuB08+ Vax/unvax v. unvax B08- (group 1,2,3)
  targetfile$bioreps <- paste(targetfile$Study_group, targetfile$Timepoint, sep = "_")
  #Unvax MamuB08+ EC v. CP (group 2)
  #targetfile$bioreps <- paste(targetfile$EC_status, targetfile$Timepoint, sep = "_")
  #Vax MamuB08+ EC v. CP (group 1)
  #targetfile$bioreps <- paste(targetfile$EC_status, targetfile$Timepoint, sep = "_")

  bioreps <- factor(targetfile$bioreps)
  biorepsdesign <- model.matrix(~ 0 + bioreps)

  # normalize
  CM2 <- DGEList(counts = CM2)
  CM2 <- calcNormFactors(CM2, method = "TMM") # TMM normalization
  message("STATUS: normalizing")
  png(file.path(norm_results, "mean_variance_norm.png"))
  Pi.CPM <- voom(counts = CM2, design=biorepsdesign,
                 normalize.method = "none", plot = T, span = 0.1, save.plot=T)
  dev.off()
  write.csv(Pi.CPM$E, file.path(norm_results, paste0("1.norm_matrix_", grodwn, ".csv")))
  message("STATUS: getting the corfit")
  corfit <- duplicateCorrelation(CM2$counts,
                                 block = factor(targetfile$Animal_ID)
  ) # account for repeated sampling of individuals

  message("STATUS: renormalizing with corfit")
  png(file.path(norm_results, "mean_variance_norm_corfit.png"))
  Pi.CPM <- voom(counts = CM2, normalize.method = "none", design=biorepsdesign,
                 correlation=corfit$consensus.correlation,
                 plot = T, span = 0.1, save.plot=T)
  dev.off()

  message("STATUS: recalculating corfit")
  corfit <- duplicateCorrelation(Pi.CPM,
                                 block = factor(targetfile$Animal_ID)
  )
  write.csv(Pi.CPM$E, file.path(norm_results, paste0("1.norm_matrix_", grodwn, "_corfit.csv")))

  sig_HGNC <- merge(rhesus2human, Pi.CPM$E,
                    by.x = "Gene.stable.ID",
                    by.y = "row.names",
                    all.X = T, all.Y = T
  )

  sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
  sig_HGNC <- avereps(sig_HGNC,
                      ID = sig_HGNC$HGNC.symbol
  )
  rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
  sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
  sig_HGNC <- as.matrix(data.frame(sig_HGNC))
  write.csv(sig_HGNC, file.path(norm_results, paste0("1.norm_matrix_HGNC_", grodwn, ".csv")), quote = FALSE)
  return(list("norm"=Pi.CPM, "corfit"=corfit))
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca=FALSE, legend="right",
                    morePCs=FALSE) {

  # Run PCA/SVD reduction
  if (isFALSE(pca)) {
    pca <- prcomp(t(exprs))
  }
  E <- get_eig(pca)
  cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
  sv <- svd(cx)


  vizualize_pca(
    file.path(results_path, paste0("svd_", base_file_name)),
    sv$u, labels[, target_columns[1]],
    labels[, target_columns[2]], figres, E, size, legend
  )
  vizualize_pca(
    file.path(results_path, paste0("pca_", base_file_name)),
    pca$x, labels[, target_columns[1]],
    labels[, target_columns[2]],
    figres, E, size, legend
  )
  vizualize_scree_plot(
    file.path(
      results_path,
      paste0("scree_", base_file_name)
    ), pca, figres
  )

  loadingscores <- as.data.frame(pca$rotation)
  is_pc1_0 <- loadingscores$PC1 > 0
  is_pc2_0 <- loadingscores$PC2 > 0

  loadingscores <- loadingscores[is_pc1_0, ]
  loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
  save_loading_scores(
    file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
    loadingscores["PC1"], figres
  )

  loadingscores <- as.data.frame(pca$rotation)
  loadingscores <- loadingscores[is_pc2_0, ]
  loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
  save_loading_scores(
    file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
    loadingscores["PC2"], figres
  )
  if (isTRUE(morePCs)) {
    is_pc3_0 <- loadingscores$PC3 > 0
    is_pc4_0 <- loadingscores$PC4 > 0
    is_pc5_0 <- loadingscores$PC5 > 0
    is_pc5_0 <- loadingscores$PC5 > 0
    is_pc6_0 <- loadingscores$PC6 > 0
    is_pc7_0 <- loadingscores$PC6 > 0

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc3_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC3)), ]
    save_loading_scores(
      file.path(results_path, paste0("loadingscores_pc3", base_file_name, ".txt")),
      loadingscores["PC3"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc4_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC4)), ]
    save_loading_scores(
      file.path(results_path, paste0("loadingscores_pc4", base_file_name, ".txt")),
      loadingscores["PC4"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc5_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC5)), ]
    save_loading_scores(
      file.path(results_path, paste0("loadingscores_pc5", base_file_name, ".txt")),
      loadingscores["PC5"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc6_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC6)), ]
    save_loading_scores(
      file.path(results_path, paste0("loadingscores_pc6", base_file_name, ".txt")),
      loadingscores["PC6"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc7_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC7)), ]
    save_loading_scores(
      file.path(results_path, paste0("loadingscores_pc7", base_file_name, ".txt")),
      loadingscores["PC7"], figres
    )

  }
  return(pca)
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, UMAP=FALSE, legend="right") {
  # Runs default paramaters of umap
  if (isFALSE(UMAP)) {
    UMAP <- umap(t(exprs))
  }
  vizualize_umap(
    file.path(results_path, paste0("umap_", base_file_name)),
    UMAP$layout, labels[, target_columns[1]],
    labels[, target_columns[2]], figres, size, legend
  )

  return(UMAP)
}

vizualize_umap <- function(plot_file, U, class1, class2, figres, size, legend) {
  # Vizualize umap reduction
  library(Polychrome)
  minx <- min(U[, 1])
  maxx <- max(U[, 1])
  miny <- min(U[, 2])
  maxy <- max(U[, 2])
  P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))

  qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
    theme_Publication() + theme(legend.title = element_blank()) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    scale_color_manual(values = as.character(P36)) +
    scale_fill_manual(values = as.character(P36)) +
    xlim(minx, maxx) + ylim(miny, maxy) +
    theme(legend.position = legend) +
    scale_shape_manual(values = seq(1, length(levels(factor(class1)))))

  ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size, legend) {
  # Vizualize PCA  results
  library(Polychrome)
  minx <- min(PCA[, 1])
  maxx <- max(PCA[, 1])
  miny <- min(PCA[, 2])
  maxy <- max(PCA[, 2])

  P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
  qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
    theme_Publication() +
    theme(legend.title = element_blank()) +
    xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
    ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
    theme(legend.position = legend) +
    scale_color_manual(values = as.character(P36)) +
    scale_fill_manual(values = as.character(P36)) +
    scale_shape_manual(values = seq(1, length(levels(factor(class1)))))

  ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

vizualize_scree_plot <- function(plot_file, PCA, figres) {
  # Vizualize principle component variation results
  scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
  png(plot_file, width = 7, height = 6, units = "in", res = figres)
  print(scree.plot)
  dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
  # Save list of genes that have a positive effect on variation of principle
  # component 1 and 2 sorted from most influential
  write.table(df, file = write_file)
}

filter_read_counts_mean <- function(cm, filter_cutoff) {
  # Filter value was calculated by:
  #### Filters by row means, test iteratively. Commonly 10 reads per gene across all samples

  A <- rowMeans(cm)
  isexpr <- A >= filter_cutoff
  cmfl <- cm[isexpr, ]
  return(cmfl)
}

rename_samples <- function(samples) {
  newsampleIDs <- c()
  for (i in samples) {
    i <- str_remove(i, "_RNA\\S+_Lib1\\S*$")
    i <- str_replace_all(i, "-", "_")
    newsampleIDs <- c(newsampleIDs, i)
  }
  return(newsampleIDs)
}

convertensembl2hgnc <- function(df) {
  df_HGNC <- merge(rhesus2human, df,
                   by.x = "Gene.stable.ID",
                   by.y = "row.names",
                   all.X = T, all.Y = T
  )

  df_HGNC <- df_HGNC[, !(names(df_HGNC) %in% c("Gene.stable.ID"))]
  df_HGNC <- avereps(df_HGNC,
                     ID = df_HGNC$HGNC.symbol
  )
  rownames(df_HGNC) <- df_HGNC[, "HGNC.symbol"]
  df_HGNC <- df_HGNC[, !(colnames(df_HGNC) %in% c("HGNC.symbol"))]
  df_HGNC <- as.matrix(data.frame(df_HGNC, check.names = FALSE))
  class(df_HGNC) <- "numeric"
  return(df_HGNC)
}

generate_design_matrix <- function(normmatrix, target) {
  rep <- factor(target$bioreps)
  sex <- factor(target$Sex)
  mm <- model.matrix(~ 0 + rep)

  rownames(mm) <- colnames(normmatrix)
  colnames(mm) <- make.names(colnames(mm))
  mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
  mm <- mm[, colSums(mm) > 0]

  excludeAll <- nonEstimable(mm)
  if (length(excludeAll) > 0) {
    message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
  }

  if ("ti" %in% excludeAll) {
    return("interactions term non estimable")
  }
  mm <- mm[, !colnames(mm) %in% excludeAll]
  if (!is.fullrank(mm)) {
    return("not full rank")
  }
  return(mm)
}

########## -- Load Data --##########

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("../count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("../target.csv", sep = ",",
                   row.names = 1, as.is = TRUE,  header = TRUE, stringsAsFactors = FALSE)

# --Rename sample names in count matrix

if (all.equal(colnames(cm), rownames(target)) != TRUE) {
  print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
  cm <- cm[, rownames(target)]
  if (all.equal(colnames(cm), rownames(target)) == TRUE) {
    print("ISSUE HAS BEEN FIXED")
  } else {
    print("ISSUE IS NOT FIXED, PLEASE LOOK AT MANUALLY")
  }
}

# --Write new renamed files to file
count_results <- "1.count_data/"
generate_folder(count_results)
write.csv(cm, file.path(count_results, "count_matrix_renamed.csv"))
write.csv(target, file.path(count_results, "target_renamed.csv"))

# --generate figure of all counts
generate_density_plot(
  cm, rownames(target), file.path(count_results, "de_intensities_raw_counts.png"), 100
)

# -- Normalize data
norm_results <- "1.norm_data/"
generate_folder(norm_results)
nd <- normalize_data(cm, target, "all")
Pi.CPM = nd$norm
corfit = nd$corfit
generate_boxplots_voom(Pi.CPM$E, target,
                       file.path(norm_results, "boxplot_vnorm_all.png"),
                       100,
                       maintitle = "Normalized count matrix",
                       ylabtitle = "voom normalized expression"
)
generate_density_plot(
  Pi.CPM$E, rownames(target), file.path(norm_results, "de_intensities_norm_counts.png"), 100
)
#-- filter out genes from each grodwn that are below mean count of 10 across samples
#Iteratively adjusted thresholds and decided that 3 was the best cutoff to get rid of the
# mean-variance hook shown in the initial mean-variance plot in 1.norm_results/ (iterative results not saved just ran in R)
cmfl_counts <- filter_read_counts_mean(cm, 3)
write.csv(cmfl_counts, file.path(count_results, "count_matrix_renamed_fl.csv"))

norm_results <- "1.norm_data_fl/"
generate_folder(norm_results)
nd2 <- normalize_data(cmfl_counts, target, "all")
Pi.CPM = nd2$norm
corfit = nd2$corfit
generate_boxplots_voom(Pi.CPM$E, target,
                       file.path(norm_results, "boxplot_vnorm_all.png"),
                       100,
                       maintitle = "Normalized count matrix",
                       ylabtitle = "voom normalized expression"
)

generate_density_plot(
  Pi.CPM$E, rownames(target), file.path(norm_results, "de_intensities_norm_counts.png"), 100
)
saveRDS(Pi.CPM, file.path(norm_results, "normobject.rds"))
saveRDS(corfit, file.path(norm_results, "corfitobject.rds"))

########## -- Feature Reduction --##########

message("STATUS: Run Feature reduction")
feature_results <- "1.feature_red"
generate_folder(feature_results)
pca <- pca_fun(
  Pi.CPM$E, target,
  feature_results, "_norm.png",
  c("Study_group", "Timepoint"), 300, 3
)

pca <- pca_fun(
  Pi.CPM$E, target,
  feature_results, "_sexnorm.png",
  c("Sex", "Study_group"), 300, 3
)

umap <- umap_fun(
  Pi.CPM$E, target,
  feature_results, "_norm.png",
  c("Timepoint", "Animal_ID"), 300, 3,
)


########## -- run lmfit --##########
message("STATUS: Running Lmfit")

#MamuB08+ Vax/unvax v. unvax B08- (group 1,2,3)
target$bioreps <- paste(target$Study_group, target$Timepoint, sep = "_")
#Unvax MamuB08+ EC v. CP (group 2)
#target$bioreps <- paste(target$EC_status, target$Timepoint, sep = "_")
#Vax MamuB08+ EC v. CP (group 1)
#target$bioreps <- paste(target$EC_status, target$Timepoint, sep = "_")
mm_all <- generate_design_matrix(Pi.CPM$E, target)
deresults_path <- "1.de"
generate_folder(deresults_path)
Pi.lmfit <- lmFit(Pi.CPM, design = mm_all, block=target$Animal_ID,
                  correlation=corfit$consensus)
saveRDS(Pi.lmfit, file.path(deresults_path, "lmfitobject.RDS"))
# Pi.lmfit <- readRDS( file.path(deresults_path, "lmfitobject.RDS"))

########## -- de analysis --##########
message("STATUS: Running DE")


#MamuB08+ Vax/unvax v. unvax B08- (group 1,2,3)
contrastsmatrix <- c(
  "(repGroup_1_D3 -repGroup_1_D0)",
  "(repGroup_1_D7 -repGroup_1_D0)",
  "(repGroup_1_D10 -repGroup_1_D0)",
  "(repGroup_1_D14 -repGroup_1_D0)",
  "(repGroup_2_D3 -repGroup_2_D0)",
  "(repGroup_2_D7 -repGroup_2_D0)",
  "(repGroup_2_D10 -repGroup_2_D0)",
  "(repGroup_2_D14 -repGroup_2_D0)",
  "(repGroup_3_D3 -repGroup_3_D0)",
  "(repGroup_3_D7 -repGroup_3_D0)",
  "(repGroup_3_D10 -repGroup_3_D0)",
  "(repGroup_3_D14 -repGroup_3_D0)"
)

  #Unvax MamuB08+ EC v. CP (group 2) and Vax MamuB08+ EC v. CP (group 1) - different sample input
  # contrastsmatrix <- c(
  # "(repEC_D3 -repEC_D0)",
  #	"(repEC_D7 -repEC_D0)",
  #	"(repEC_D10 -repEC_D0)",
  #	"(repEC_D14 -repEC_D0)",
  #	"(repNonEC_D3 -repNonEC_D0)",
  #	"(repNonEC_D7 -repNonEC_D0)",
  #	"(repNonEC_D10 -repNonEC_D0)",
  #	"(repNonEC_D14 -repNonEC_D0)"
  # )


contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm_all)
contrast <- contrasts.fit(Pi.lmfit, contrasts = contr)
Pi.contrasts <- eBayes(contrast, robust = TRUE, trend = TRUE)
results <- decideTests(Pi.contrasts, lfc = 0.58, method = "separate", adjust.method = "BH", p.value = 0.05)
write.csv(Pi.contrasts$coefficients, file = file.path(deresults_path, "coefficients.csv"), quote = F)
write.csv(Pi.contrasts$t, file = file.path(deresults_path, "t_stats.csv"), quote = F)
write.csv(Pi.contrasts$p.value, file = file.path(deresults_path, "p_value.csv"), quote = F)
for (i in 1:ncol(Pi.contrasts$p.value)) Pi.contrasts$p.value[, i] <- p.adjust(Pi.contrasts$p.value[, i], method = "BH")
write.csv(Pi.contrasts$p.value, file = file.path(deresults_path, "p_value_adj.csv"), quote = F)
dataMatrixde <- Pi.contrasts$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)
Pi.contrasts$genes <- data.frame(ID_REF=rownames(Pi.contrasts))
write.csv(ExpressMatrixde, file = file.path(deresults_path, "expression_matrix_de.csv"), quote = F)
write.csv(results, file = file.path(deresults_path, "results_de.csv"), quote = F)
write.csv(dataMatrixde, file = file.path(deresults_path, "full_expression_matrix_de.csv"), quote = F)

#Formatted input for bargraph fig

ggplot(DE_bargraph)
  + geom_bar(aes(as.factor(Day),Count, fill=Variable), position="stack", stat = "identity", colour="black")
  + theme_bw()
  + scale_fill_manual(values = c("#0571b0","#ca0020"))
  + xlab("Day")
  + ylab("Number of DE genes")
  + facet_wrap(~Group, ncol=1)


########## -- dde analysis --##########

message("STATUS: Running DDE")

#Same contrasts for all ddes, different sample input by group and controller status
#See metadata in target file
contrastsmatrix <- c(
  "(repGroup_1_D3 - repGroup_1_D0) - (repGroup_2_D3 - repGroup_2_D0)",
  "(repGroup_1_D7 - repGroup_1_D0) - (repGroup_2_D7 - repGroup_2_D0)",
  "(repGroup_1_D10 - repGroup_1_D0) - (repGroup_2_D10 - repGroup_2_D0)",
  "(repGroup_1_D14 - repGroup_1_D0) - (repGroup_2_D14 - repGroup_2_D0)"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm_all)
contrast <- contrasts.fit(Pi.lmfit, contrasts = contr)
Pi.contrasts <- eBayes(contrast, robust = TRUE, trend = TRUE)
results <- decideTests(Pi.contrasts, lfc = 0.58, method = "separate", adjust.method = "BH", p.value = 0.05)
write.csv(Pi.contrasts$coefficients, file = file.path(deresults_path, "coefficients.csv"), quote = F)
write.csv(Pi.contrasts$t, file = file.path(deresults_path, "t_stats.csv"), quote = F)
write.csv(Pi.contrasts$p.value, file = file.path(deresults_path, "p_value.csv"), quote = F)
for (i in 1:ncol(Pi.contrasts$p.value)) Pi.contrasts$p.value[, i] <- p.adjust(Pi.contrasts$p.value[, i], method = "BH")
write.csv(Pi.contrasts$p.value, file = file.path(deresults_path, "p_value_adj.csv"), quote = F)
dataMatrixde <- Pi.contrasts$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)
Pi.contrasts$genes <- data.frame(ID_REF=rownames(Pi.contrasts))
write.csv(ExpressMatrixde, file = file.path(deresults_path, "expression_matrix_de.csv"), quote = F)
write.csv(results, file = file.path(deresults_path, "results_de.csv"), quote = F)
write.csv(dataMatrixde, file = file.path(deresults_path, "full_expression_matrix_de.csv"), quote = F)

###For significant genes, pull and format ExpressMatrixde and p_value_adj
###For all genes or less strict cutoffs, pull and format dataMatrixde and p_value <- this is in supp data

HM <- Heatmap(as.matrix(heatmatrix), clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2", cluster_columns = FALSE, cluster_rows = TRUE, row_names_gp = gpar(fontsize = 4), name="Log2FC", col=mycols)
HM <- draw(HM)  #Sanity check or print

#Pull out genes in module/clusters for GO Enrichment analysis
r.dend <- row_dend(HM)
rcl.list <- row_order(HM)
lapply(rcl.list, function(x) length(x))
for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(as.matrix(heatmatrix[row_order(HM)[[i]],]))))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(as.matrix(heatmatrix[row_order(HM)[[i]],]))))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
