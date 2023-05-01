#' detect_bias()
#'
#' @param path = "/logfcs"
#' @param locusInfo = TRUE
#' @param savePlot = TRUE
#'
#' @return
#' @export
#'
#' @examples
#' detect_bias()
detect_bias <- function(path, locusInfo = TRUE, savePlot = TRUE){
  myfiles <- lapply(list.files(path = path, pattern = "*.csv", full.names = TRUE), read.delim, sep = ",")
  joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
  filenames <- list.files(path = path, pattern = "*.csv") %>%
    gsub(pattern = ".csv", replacement = "")
  if (locusInfo == TRUE){
    info <- joined[,c(1:3)]
    colnames(info) <- c("locus_tag", "gene", "function")
    dat <- joined[,-c(1:3)]
  } else {
    info <- as.data.frame(joined[,1])
    dat <- joined[,-c(1)]
    colnames(dat) <- "locus_tag"
  }
  dat.2 <- dat %>% select(contains(c("logFC", 'pvalue')))
  colnames(dat.2) <- c(gsub("$", "_logFC", filenames), gsub("$", "_pvalue", filenames))

  cols <- c("logFC", "pvalue", "locus_tag", "sig", "ob", "cond")
  dat.plot <- data.frame()
  for (i in 1:length(filenames)) {
    dat.input <- dat.2[,c(i, i+length(filenames))]
    dat.input$locus_tag <- info$locus_tag
    dat.input$sig <- ifelse(dat.input[,2]<0.05, "yes", "no")
    dat.input$ob <- 1:nrow(dat.input)
    dat.input$cond <- filenames[i]
    colnames(dat.input) <- cols
    dat.plot <- rbind(dat.plot, dat.input)
  }

  plottheme <- theme_bw() +
    theme(text = element_text(size = 20))

  p <- ggplot(dat.plot, aes(x = ob, y = logFC, col = sig)) +
    geom_point(size = 0.1) +
    facet_wrap(~cond) +
    labs(y = "Log2 Fold Change", x = "Locus", col = "Significant?") +
    scale_color_manual(values = c("yes" = "red", "no" = "black")) +
    scale_x_continuous(breaks = c(0, (round(nrow(info), digits = -3))/2, round(nrow(info), digits = -3))) +
    plottheme

  if (savePlot == TRUE){
    png(filename = paste0(path, "detect_bias.png"), res = 300, height = 2500, width = 3500)
    suppressWarnings(print(p))
    dev.off()
    paste0("Your diagnostics file has been saved to ", path, "detect_bias.png")
  } else {
    suppressWarnings(plot(p))
  }
}
