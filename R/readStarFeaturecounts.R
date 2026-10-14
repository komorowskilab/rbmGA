#' @title readStarFeaturecounts
#' @author Girish Pulinkala
#' @description Make a Decision table from star read counts,The last column of decision should be missing at this point. For experts, this is an information system.
#' @param folderpath  A string value pointing to the folder containing read counts for each sample individually for .e.g. 531e82fc-c135-4700-9887-3ea845ce8b2c.rna_seq.augmented_star_gene_counts.tsv
#' @param makedt A boolean value that specifies whether to generate a decision table or a dataframe for differential expression.
#' @returns A dataframe for ML or Differential expression analysis
#' @export
#' @import data.table
#' @import purrr
#' @import dplyr
#' @import tidyr
#' @importFrom utils read.delim
#' @examples
#' # path <- c("../data/data_all_old/")
#' # df<- readFeaturecounts(path,makedt=TRUE)

readStarFeaturecounts<- function(folderpath,makedt=TRUE,  strand="unstranded",   col.names = c("gene_id", "gene_name", "gene_type",
                                                                                             "unstranded", "stranded_first", "stranded_second",
                                                                                             "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"), dtgeneid='gene_name'){

  strand <- strand   # or "stranded_first", "stranded_second", "tpm_unstranded", "reverse", ...

  files <- list.files(folderpath, pattern = "(?i)\\.tsv(%0D)?$", full.names = TRUE)
  files <- set_names(files, sub("(?i)\\.tsv(%0D)?$", "", basename(files), perl = TRUE))

  star_long <- map_dfr(
    files,
    ~ read.delim(.x, header = TRUE, sep = "\t", comment.char = "#",
                 col.names = col.names,
                 stringsAsFactors = FALSE) %>%
      dplyr::filter(!grepl("^N_", trimws(gene_id))),
    .id = "sample"
  )

  df <- star_long %>%
    dplyr::select(all_of(dtgeneid), sample, count = all_of(strand)) %>%
    tidyr::pivot_wider(names_from = sample, values_from = count, values_fn = sum)


  df<- df %>% tibble::column_to_rownames(.,dtgeneid)

  if(makedt==TRUE){
    df <- as.data.frame(t(df))
  }

  return(as.data.frame(df))
}

