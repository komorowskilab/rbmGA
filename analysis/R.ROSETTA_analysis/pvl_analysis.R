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

pvl_data<- fread('../../data/external_validation/pc531_36514 viral load data.csv')
pvl_data$PID<- as.factor(pvl_data$PID)
pvl_data$`F5=Rh36514 (P)`<-gsub(pattern = ',','',pvl_data$`F5=Rh36514 (P)`)
pvl_data$`F5=Rh36514 (P)`<-as.numeric(pvl_data$`F5=Rh36514 (P)`)

pvl_data %>%
  filter(!is.na(`F5=Rh36514 (P)`)) %>%
  ggplot(aes(x = factor(PID), y = `F5=Rh36514 (P)`)) +
  geom_line(aes(group = 1),size = 1,color = "red") +
  geom_point(size = 3, color='red') +
  labs(x = "PID (Post Infection Day)", y='Plasma Viral Load (PVL) after 5xChallenge')+ theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


