#' Uncorrected ciprofloxacin transposon insertion sequencing dataset affected by chromosomal location bias
#'
#' @format
#' A data frame with 3,920 rows and 7 columns:
#' \describe{
#'   \item{locus_tag}{Escherichia coli BW25113 locus tag}
#'   \item{gene_name}{Gene name of the associated locus_tag}
#'   \item{function.}{Function of the associated locus_tag}
#'   \item{logFC}{Log2 fold change calculated with the Bio::TraDIS toolkit, without normalisation}
#'   \item{logCPM}{Log2 counts per million calculated with the Bio::TraDIS toolkit, without normalisation}
#'   \item{PValue}{Pvalue (unadjusted) calculated with the Bio::TraDIS toolkit, without normalisation}
#'   \item{q.value}{Adjusted pvalue calculated with the Bio::TraDIS toolkit, without normalisation}
#' }
#' @source <paper link will go here>
"cip_uncorrected"
