\name{structure_rc}
\alias{structure_rc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine read counts from multiple files listed in a folder
}
\description{
This function takes all .csv files in a specified directory 
(default: current working directory) and merges them by column 'locus_tag.' 
It then extracts the column called 'read_count' and puts them together in a 
data frame. 
}
\usage{
structure_rc(csvpath = getwd(), getLocusInfo = TRUE, suffix = ".tradis_gene_insert_sites.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csvpath}{
Path to folder containing .csv files containing read counts. All .csv files will be imported. \cr
Default: current working directory.
}
  \item{getLocusInfo}{
TRUE or FALSE. If TRUE, locus information is saved as 'locusInfo.tsv' for import during normalisation. If FALSE, normalisation won't be able to run for normalisation and only locus tags will be outputted. \cr
Default: TRUE.
  }
    \item{suffix}{
Read count files are expected to be traDIS outputs, so expected names are "condition_rep.tradis_gene_insert_sites.csv". If your file type is different, change to everything after the replicate information to rename columns efficiently.\cr
Default: .tradis_gene_insert_sites.csv
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Reference goes here:
}
\author{
Geraldine J Sullivan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# If current working directory is folder with all .csv files:
readcounts <- structure_rc()

# If wanting to specify path:
readcounts <- structure_rc(csvpath = "~/path/to/files")

# If not wanting locus information:
readcounts <- structure_rc(csvpath = "~/path/to/files", getLocusInfo = FALSE)

