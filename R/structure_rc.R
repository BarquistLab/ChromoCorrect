#' @title Structure read counts
#'
#' @description Structure read counts for normalisation with normalise_bias()
#'
#' @author Geraldine J Sullivan
#'
#' @param csvpath Path to folder containing .csv files containing read counts. All .csv files will be imported. Default: current working directory.
#' @param getLocusInfo TRUE or FALSE. If TRUE, locus information is saved as 'locusInfo.tsv' for import during normalisation. If FALSE, normalisation won't be able to run for normalisation and only locus tags will be outputted. Default: TRUE.
#' @param suffix Read count files are expected to be traDIS outputs, so expected names are "condition_rep.tradis_gene_insert_sites.csv". If your file type is different, change to everything after the replicate information to rename columns efficiently. Default: .tradis_gene_insert_sites.csv
#'
#' @references Will go here later
#'
structure_rc <- function(csvpath = "/readcounts", getLocusInfo = TRUE, suffix = ".tradis_gene_insert_sites.csv") {
  myfiles <- lapply(list.files(path = csvpath, pattern = "*.csv",
                               full.names = TRUE), read.delim)
  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
  filenames <- list.files(pattern = "*sites.csv", path = csvpath) %>%
    gsub(pattern = suffix, replacement = "")
  rc <- joined %>% select(contains(c("locus_tag", "read_count")))
  colnames(rc)[2:ncol(rc)] <- filenames
  if (getLocusInfo == TRUE){
    locusInfo <- joined[,c(1,2,11)]
    colnames(locusInfo) <- c("locus_tag", "gene_name", "Function")
    write.table(locusInfo, file = paste0(csvpath, "/locusInfo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("LocusInfo written to", csvpath)
  }
  return(rc)
}

