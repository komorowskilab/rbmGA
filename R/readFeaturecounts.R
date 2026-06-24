#' @title readFeaturecounts
#' @author Girish Pulinkala
#' @description Make a Decision table from htseq read counts,The last column of decision should be missing at this point. For experts, this is an information system.
#' @param folderpath  A string value pointing to the folder containing read counts for each sample individually for .e.g. SJACT001_D.RNA-Seq.feature-counts.txt
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

readFeaturecounts<- function(folderpath,makedt=TRUE){

      df <- folderpath %>% list.files(pattern = "(?i)\\.txt|.txt%0D$", full.names = TRUE) %>%
      purrr::set_names(tools::file_path_sans_ext(basename(.))) %>%
      purrr::map_dfr(utils::read.delim, col.names = c('V1', "V2"), .id = "colname") %>%
      pivot_wider(names_from = colname, values_from = V2)

      df<- df %>% tibble::column_to_rownames(.,'V1')

      if(makedt==TRUE){
        df <- as.data.frame(t(df))
      }

      return(as.data.frame(df))
}


