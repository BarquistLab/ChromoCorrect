#' @title detect_bias
#'
#' @description Generate a png of locus by fold change scatterplot for each output file in a directory to determine whether any are affected by chromosomal location bias.
#'
#' @param path = "/logFCs"
#' @param locusInfo = TRUE
#' @param savePlot = TRUE
#'
#' @export
#'
detect_bias <- function(path = "/logfcs", locusInfo = TRUE, savePlot = TRUE){
  myfiles <- lapply(list.files(path = path, pattern = "*.csv", full.names = TRUE), utils::read.delim)
  joined <- myfiles %>% purrr::reduce(dplyr::full_join, by = "locus_tag")
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
  dat.2 <- dat %>% dplyr::select(dplyr::contains(c("logFC", 'pvalue')))
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
    theme(text = ggplot2::element_text(size = 20))

  p <- ggplot2::ggplot(dat.plot, aes(x = ob, y = logFC, col = sig)) +
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

  res <- data.frame()
  for (h in 1:length(myfiles)){
    sum.dat <- myfiles[[h]]
    sum.dat <- sum.dat[!is.na(sum.dat$logFC),]
    length <- ceiling(nrow(sum.dat)/5)
    sum.dat2 <- split(sum.dat, rep(1:ceiling(nrow(sum.dat)/length), each=length, length.out=nrow(sum.dat)))

    summary <- data.frame()
    for (i in 1:length(sum.dat2)){
      split.dat <- sum.dat2[[i]]
      split.dat.2 <- cut(split.dat$logFC, quantile(split.dat$logFC, c(0, 0.2, 0.8, 1)), include.lowest = TRUE, lab = c("lo", "mid", "hi"))
      split.dat.2 <- split(split.dat$logFC, split.dat.2)
      split.dat.2 <- as.data.frame(split.dat.2$mid)
      split.dat.2$ob <- 1:nrow(split.dat.2)
      model <- lm(split.dat.2$`split.dat.2$mid`~split.dat.2$ob)
      summary[i,1] <- summary(model)$coefficients[2,4]
      summary[i,2] <- mean(split.dat.2$`split.dat.2$mid`)
      summary[i,3] <- summary(model)$coefficients[1,1]
      summary[i,4] <- max(split.dat.2$`split.dat.2$mid`)-min(split.dat.2$`split.dat.2$mid`)
    }
    summary$sample <- filenames[h]
    res <- rbind(res, summary)
  }
  if (any(res$V1<0.1) & any(abs(res$V3)>0.05)) {
    sub.dat <- res[res$V1<0.1,]
    sub.dat <- sub.dat[sub.dat$V3>0.05,]
    culprits <- unique(sub.dat$sample)
  }
  cat(sep = "", "It appears the following samples have chromosomal location bias: ", paste0(culprits, collapse = ", "), ".\nPlease consider proceeding to correction with these samples.")
}
