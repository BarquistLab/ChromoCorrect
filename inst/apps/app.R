library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyalert)
library(shinydashboardPlus)
library(htmltools)
library(htmlwidgets)
library(dplyr)
library(shinyjs)
library(ggplot2)
library(patchwork)
library(edgeR)
library(FBN)
library(shinyalert)

wdir <- getwd()

shinyApp(
  ui = dashboardPage(
    options = list(sidebarExpandOnHover = TRUE),
    header = dashboardHeader(title = "ChromoCorrect", titleWidth = 300),

    sidebar = dashboardSidebar(width = 300, minified = F, collapsed = F,
                               h4("Upload files here"),
                               uiOutput("mytab1"), uiOutput("mytab1.1"),
                               uiOutput("mytab2"), uiOutput("mytab2.1")
    ),
    body = dashboardBody(
      tags$head(tags$script('
                                var dimension = [0, 0];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerWidth-300;
                                    dimension[1] = window.innerHeight-300;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerWidth-300;
                                    dimension[1] = window.innerHeight-300;
                                    setTimeout(function() {
                                      Shiny.onInputChange("dimension", dimension)
                                    }, 500);
                                });
                            ')),
      h3('Detecting and correcting chromosomal bias'),
      tabsetPanel(id = "tabs",
                  tabPanel("Detecting",
                           br(),
                           textOutput("dimension_display"),
                           h4("Upload your Bio::TraDIS output files to determine whether chromosomal bias is affecting your data"),
                           p("If the overall trend of your fold changes does not match the red line, your data needs normalising."),
                           box(title = "Locus by fold change scatterplot",
                               status = "primary",
                               width = 8, height = 6,
                               imageOutput("detec_fc", height = "100%", width = "100%")),
                           box(
                             width = 4,
                             title = "Decision", status = "warning",
                             h4(htmlOutput(outputId = "detec_text"))
                           )),

                  tabPanel("Correcting",
                           br(),
                           h4("Upload your Bio::TraDIS read files to correct the chromosomal bias affecting your data"),
                           p("This requires two control files and two condition files, or one file of read counts containing all conditions of interest."),
                           box(title = "Before and after normalisation",
                               status = "primary",
                               width = 12, height = 6,
                               imageOutput("corrected_plot", height = "100%", width = "100%")),
                           box(title = "Normalised data",
                               uiOutput("download"), br(),
                               status = "primary",
                               width = 12, height = 6,
                               DT::dataTableOutput("normdata"))
                  ))
    )
  ),
  server = function(input, output){

    output$test <- renderText({
      return(paste("<style='text-indent:1em;'>If the overall trend of your fold changes does not match the red line, your data needs normalising."))
    })

    output$mytab1 <- renderUI({
      tagList(
        conditionalPanel(condition = 'input.tabs=="Detecting"',
                         fileInput("uploadfc", "Upload your traDIS output file(s) here", buttonLabel = "Browse...", multiple = TRUE),
        ))
    })

    output$mytab1.1 <- renderUI({
      tagList(
        conditionalPanel(condition = 'input.tabs=="Detecting"',
                         selectizeInput("datasetsnorm", "Select dataset for visualising:",
                                        choices = gsub(".csv", "", input$uploadfc$name))
        ))
    })

    detecplot <- reactive({
      req(input$uploadfc)
      num <- grep(value = FALSE, pattern = input$datasetsnorm, x = input$uploadfc$name)
      data <- read.csv(input$uploadfc[[num, "datapath"]])
      data$obs <- 1:nrow(data)
      data$`Significance (0.05)` <- ifelse(data$q.value<0.05, "Significant", "Not significant")
      ggplot(data, aes(x = obs, y = logFC, col = `Significance (0.05)`)) +
        geom_point(cex = 0.5) +
        geom_hline(yintercept = 0, col = "red") +
        theme_classic() +
        scale_color_manual(values = c("Significant" = "red", "Not significant" = "black")) +
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size = 16)) +
        labs(x = "Locus", y = "Log2 Fold Change", title = paste("Fold change by locus scatterplot - ",
                                                                gsub(pattern = ".csv", "", input$uploadfc$name[[num]])))
    })

    output$detec_fc <- renderImage({
      req(input$uploadfc)
      outfile <- tempfile(fileext = ".png")
      png(outfile,
          width = 0.6*input$dimension[1]*8,
          height = 400*8,
          res = 72*8)
      print(detecplot())
      dev.off()

      list(src = outfile,
           contentType = 'image/png',
           width = 0.6*input$dimension[1],
           height = 400,
           alt = "This is alternate text")
    }, deleteFile = TRUE)

    output$detec_text <- renderText({
      req(input$uploadfc)
      num <- grep(value = FALSE, pattern = input$datasetsnorm, x = input$uploadfc$name)
      data <- read.csv(input$uploadfc[[num, "datapath"]])
      data$obs <- 1:nrow(data)
      length <- ceiling(nrow(data)/5)
      datcut <- split(data, rep(1:ceiling(nrow(data)/length), each=length, length.out=nrow(data)))
      summary <- data.frame()
      for (i in 1:length(datcut)){
        dat <- datcut[[i]]
        dat <- dat[complete.cases(dat$logFC),]
        dat.1 <- cut(dat$logFC, quantile(dat$logFC, c(0, 0.2, 0.8, 1)), include.lowest = TRUE, lab = c("lo", "mid", "hi"))
        dat.2 <- split(dat$logFC, dat.1)
        dat.2 <- as.data.frame(dat.2$mid)
        dat.2$ob <- 1:nrow(dat.2)
        model <- lm(dat.2$`dat.2$mid`~dat.2$ob)
        summary[i,1] <- summary(model)$coefficients[2,4]
        summary[i,2] <- median(dat.2$`dat.2$mid`)
        summary[i,3] <- mean(dat.2$`dat.2$mid`)
      }
      if (any(summary$V1<0.1) & (any(abs(summary$V2)>0.2) | any(abs(summary$V3)>0.2))){
        return(paste("<span style=\"color:red\">The trend line does not appear to be equal to 0.<br>Please consider proceeding to correction.</span>"))
      } else {
        return(paste("<span style=\"color:green\">The trend line appears to be approximately equal to 0.<br>Your data does not need correction.</span>"))
      }
    })

    output$mytab2 <- renderUI({
      tagList(
        conditionalPanel(condition = 'input.tabs=="Correcting"',
                         fileInput("uploadrc", "Upload your traDIS read files here", multiple = TRUE, accept = c(".csv", ".tsv")),
                         fileInput("rcfile", "OR upload your read count table here", multiple = FALSE)
        ))
    })

    output$mytab2.1 <- renderUI({
      tagList(
        conditionalPanel(condition = 'input.tabs=="Correcting"',
                         p("Your first column should be 'locus_tag', with read counts for replicates of a control and a condition in the following columns. Replicate column names should end with '_1', '_2' etc"),
                         h4("Choose condition for control here"),
                         if (!is.null(input$uploadrc)) {
                           selectizeInput("controlrc", "Select which condition is your control:",
                                          choices = c("Select one here", unique(gsub("_[0-9].tradis.gene.insert.sites.csv", "", input$uploadrc$name))))
                         } else if (!is.null(input$rcfile)) {
                           rc <- read.delim(input$rcfile$datapath)
                           selectizeInput("controlrc", "Select which condition is your control:",
                                          choices = c("Select one here", unique(gsub("_[0-9].*", "", colnames(rc)[2:ncol(rc)]))))
                         } else {
                           selectizeInput("controlrc", "Select which condition is your control:",
                                          choices = c("Select one here"))
                         },
                         uiOutput("button"),
                         hr(), h4("Optional"),
                         numericInput("minrc", "Minimum read count cutoff", value = 10),
                         fileInput("locusinfo", "Upload tab separate locus information", multiple = FALSE),
                         p("For example, a tab separated file with locus_tag, gene_name, function for extra information in the outputs.")
        ))
    })

    output$button <- renderUI({
      if (input$controlrc == "Select one here"){
        NULL
      } else if (!is.null(input$uploadrc$datapath)) {
        if (length(input$uploadrc$datapath)<4){
          shinyalert(text = "Less than 4 files detected. Please provide at least two replicates per control/condition.",
                     type = "error")
          NULL
        } else if (length(unique(gsub("_[0-9].tradis.gene.insert.sites.csv", "", input$uploadrc$name)))>2){
          shinyalert(text = "More than 2 conditions detected. Please provide at least two repliates for one control and one condition.",
                     type = "error")
          NULL
        } else {
          actionButton("run", "Start normalisation", class = "btn-primary btn-lg")
        }
      } else if (!is.null(input$rcfile)){
        rc <- read.delim(input$rcfile$datapath)
        uniq <- unique(gsub("_[0-9]$", "", colnames(rc)))
        if (length(uniq)<3) {
          shinyalert(text = "Not enough conditions detected. Please provide at least two repliates for a control and a condition.",
                     type = "error")
        } else if (length(uniq)>3) {
          shinyalert(text = "More than 2 conditions detected. Please provide at least two repliates for one control and one condition, or check your column names.",
                     type = "error")
        } else {
          actionButton("run", "Start normalisation", class = "btn-primary btn-lg")
        }
      }
    })

    readcounts <- eventReactive(input$run, {
      if (!is.null(input$uploadrc)){
        myfiles <- purrr::map(input$uploadrc$datapath, read.delim) %>%
          purrr::set_names(input$uploadrc$name)
        joined <- myfiles %>% purrr::reduce(full_join, by = "locus_tag")
        filenames <- input$uploadrc$name %>%
          gsub(pattern = ".tradis_gene_insert_sites.csv", replacement = "")
        rc <- joined %>% select(contains(c("locus_tag", "read_count")))
        colnames(rc)[2:ncol(rc)] <- filenames
      } else if (!is.null(input$rcfile)) {
        rc <- read.delim(input$rcfile$datapath)
      }
      rc
    })

    done <- FALSE

    correction <- reactive({
      window_size = 500
      while(TRUE){
        rc <- readcounts()
        rc <- cbind("locus_tag" = rc[,1], "ob" = 1:nrow(rc), rc[,2:ncol(rc)])
        norm_counts <- data.frame(row.names = rc$locus_tag)
        for (i in 3:ncol(rc)){
          calc <- rc[,c(1, 2, i)]
          calc[(nrow(calc)+1):(nrow(calc)+1000),] <- calc[1:1000,]
          calc$keep <- 1:nrow(calc)
          calc$pred <- FBN::medianFilter(inputData = calc[,3], windowSize = window_size)
          calc <- calc[!calc$keep > nrow(rc),]
          calc$ratio <- calc$pred/mean(calc$pred)
          calc$norm <- as.integer(round(calc[,3]/calc$ratio))
          norm_counts[,i-2] <- calc$norm
        }

        offset <- (log(norm_counts + 0.01) - log(rc[,3:(ncol(rc))] + 0.01))
        eff.lib <- calcNormFactors(norm_counts) * colSums(norm_counts)
        offset <- sweep(offset, 2, log(eff.lib), "-")
        colnames(offset) <- colnames(norm_counts) <- colnames(rc)[3:ncol(rc)]

        # norm_counts <- norm_counts[apply(apply(rc[3:ncol(rc)], 1, ">", 10), 2, any),]
        # offset <- offset[apply(apply(rc[3:ncol(rc)], 1, ">", 10), 2, any),]
        # rc <- rc[apply(apply(rc[3:ncol(rc)], 1, ">", 10), 2, any),]
        norm_counts <- norm_counts[apply(apply(rc[3:ncol(rc)], 1, ">", input$minrc), 2, any),]
        offset <- offset[apply(apply(rc[3:ncol(rc)], 1, ">", input$minrc), 2, any),]
        rc <- rc[apply(apply(rc[3:ncol(rc)], 1, ">", input$minrc), 2, any),]

        rownames(offset) <- rownames(rc) <- rc$locus_tag
        rc <- rc[,-c(1:2)]

        # Step 3 - edgeR (differential expression)
        group <- gsub("_[0-9]", replacement = "", x = colnames(rc))
        conds_edgeR <- as.factor(unique(group))
        #ctrl <- "MH"
        ctrl <- as.character(input$controlrc)
        conds_edgeR <- relevel(conds_edgeR, ref=ctrl)
        condition <- as.character(conds_edgeR[!conds_edgeR %in% ctrl])
        design <- model.matrix(~0+group)
        contrast <- makeContrasts(contrasts = paste0("group", condition, " - group", ctrl), levels = design)
        y <- DGEList(counts=rc, group=group, genes=rownames(rc))
        y <- scaleOffset(y, -as.matrix(offset))
        y <- estimateGLMCommonDisp(y, design)
        y <- estimateGLMTagwiseDisp(y, design)
        fit <- glmFit(y, design, robust=TRUE)
        lrt <- glmLRT(fit, contrast=contrast)
        tags_after <- lrt$table

        length <- ceiling(nrow(tags_after)/5)
        tagplot <- split(tags_after, rep(1:ceiling(nrow(tags_after)/length), each=length, length.out=nrow(tags_after)))

        summary <- data.frame()
        for (i in 1:length(tagplot)){
          dat <- tagplot[[i]]
          dat.2 <- cut(dat$logFC, quantile(dat$logFC, c(0, 0.2, 0.8, 1)), include.lowest = TRUE, lab = c("lo", "mid", "hi"))
          dat.2 <- split(dat$logFC, dat.2)
          dat.2 <- as.data.frame(dat.2$mid)
          dat.2$ob <- 1:nrow(dat.2)
          model <- lm(dat.2$`dat.2$mid`~dat.2$ob)
          summary[i,1] <- summary(model)$coefficients[2,4]
          summary[i,2] <- mean(dat.2$`dat.2$mid`)
          summary[i,3] <- summary(model)$coefficients[1,1]
          summary[i,4] <- max(dat.2$`dat.2$mid`)-min(dat.2$`dat.2$mid`)
        }
        #if (any(summary$V1<0.1) & (any(abs(summary$V2)>0.05) | any(abs(summary$V3)>0.05))){
        if (any(summary$V1<0.1) & any(abs(summary$V3)>0.05)) {
          if (window_size == 200){
            showNotification(paste("Window size of 200 is minimum to retain biological significance with operons. Stopping here."))
            break
          } else {
            showNotification(paste("Window size of ", window_size, " not correct, recomputing"))
            window_size <- window_size-100
            rm(norm_counts, offset, tags_after, tagplot, summary)
          }
        } else {
          showNotification(paste("Window size of ", window_size, " correct, finishing up"))
          break
        }
      }
      d <- DGEList(counts = rc, group=group)
      plotMDS.DGEList(d, labels=group)
      d <- calcNormFactors(d)
      d <- estimateCommonDisp(d)
      d <- estimateTagwiseDisp(d)
      de.tgw <- exactTest(d,pair=c(ctrl, condition))
      #de.tgw <- exactTest(d,pair=c("MH", "Cip"))
      tags_before <- data.frame(de.tgw$table)
      tags_before$`Significance (0.05)` <- ifelse(tags_before$PValue < 0.05, "Significant", "Not significant")
      tags_before <<- tags_before
      tags_after$ob <- 1:nrow(tags_after)
      tags_after$`Significance (0.05)` <- ifelse(tags_after$PValue < 0.05, "Significant", "Not significant")
      cond <<- condition
      tags_before <<- tags_before
      tags_after <<- tags_after
      window_size <<- as.numeric(window_size)
    })

    corplot <- reactive({
      req(input$run)
      correction()
      before <- ggplot(tags_before, aes(x = 1:nrow(tags_before), y = logFC, col = `Significance (0.05)`)) +
        geom_point(size = 0.5) +
        theme_classic() +
        theme(text = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 18)) +
        labs(title = paste0("Uncorrected fold change plot by locus - ", cond),
             x = "Locus", y = "Log2 fold change") +
        scale_color_manual(values = c("Not significant" = "black", "Significant" = "red"), guide = "none")
      after <- ggplot(tags_after, aes(x = 1:nrow(tags_after), y = logFC, col = `Significance (0.05)`)) +
        geom_point(size = 0.5) +
        theme_classic() +
        theme(text = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 18),
              plot.subtitle = element_text(hjust = 0.5, size = 14)) +
        labs(title = paste0("Corrected fold change plot by locus - ", cond),
             x = "Locus", y = NULL, col = "Significant (0.05)", subtitle = paste("Sliding median window size of", as.numeric(window_size))) +
        scale_color_manual(values = c("Not significant" = "black", "Significant" = "red"))
      done <<- TRUE
      before+after
    })

    output$corrected_plot <- renderImage({
      req(input$run)
      outfile <- tempfile(fileext = ".png")
      png(outfile,
          width = 0.95*input$dimension[1]*8,
          height = 500*8,
          res = 72*8)
      print(corplot())
      dev.off()

      list(src = outfile,
           contentType = 'image/png',
           width = 0.95*input$dimension[1],
           height = 500,
           alt = "This is alternate text")
    }, deleteFile = TRUE)

    output$normdata <- DT::renderDataTable({
      req(input$run)
      output$download <- renderUI({
        actionButton("download_attempt", label = "download csv", class = "btn-secondary")
      })
      corplot()
      if (done == TRUE) {
        tags_after <- cbind("locus_tag" = rownames(tags_after), tags_after[,c(1:(ncol(tags_after)-2))])
        rownames(tags_after) <- 1:nrow(tags_after)
        colnames(tags_after)[4:5] <- c("Pvalue", "q.value")
        table_out <<- tags_after
        output$download <- renderUI(actionButton("download_attempt", "Download csv"))
        if (!is.null(input$locusinfo)){
          locusinfo <- read.delim(input$locusinfo$datapath)
          table_out <<- merge(locusinfo, table_out, by = "locus_tag", all.x = TRUE)
        }
        table_out <- table_out[order(table_out$q.value, decreasing = FALSE),]
        DT::datatable(table_out, rownames = FALSE, options = list(pageLength = 15)) %>% DT::formatRound(columns = c((ncol(table_out)-3):(ncol(table_out))), digits = c(2,2,4,4))
      }
    })

    observeEvent(input$download_attempt, {
      write.table(table_out,file=paste0(wdir, "/", cond, "_ChromoCorrect.csv"),
                  append=FALSE, quote=TRUE, sep=",", row.names=FALSE)
      shinyalert(title = "Success",
                 text = paste0(cond, "_ChromoCorrect.csv has been saved to ", wdir))
    })
  }
)
