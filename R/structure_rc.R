structure_rc <- function(csvpath = getwd(), getLocusInfo = TRUE, suffix = ".tradis_gene_insert_sites.csv") {
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

