\name{normalise_CB}
\alias{normalise_CB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalise chromosomal location bias in read counts
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Readcounts from structure_rc(). \cr
Required.
}
  \item{control}{
Name of the condition that is to be used for the control during comparisons. No replicate information desired. \cr
Required.
}
\item{windowSize}{
Specify the window size for the sliding median normalisation. Either "auto" or a number. Minimum 200, maximum 1000. \cr
Optional. Default: auto.
}
\item{minrc}{
Minimum read count for filtering out low values.  \cr
Optional. Default: 10.
}
\item{locusInfo}{
TRUE or FALSE. Whether locus information should be added to locus tags. Only possible if structure_rc() set getLocusInfo to TRUE. Path to file must be specified in path variable. \cr
Optional. Default: TRUE.
}
\item{path}{
Path to folder containing locusInfo.tsv, and path where plots will be saved. \cr
Required if locusInfo = TRUE. Default: getwd().
}
\item{writePlots}{
TRUE or FALSE. Write out locus by fold change scatterplots to observe peaks and troughs in data. Plots are written to path specified in path variable. \cr
Optional. Default: TRUE.
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
# If column names for control are "MH_1" and "MH_2":
normalise_CB(x, control = "MH", windowSize = "auto", minrc = 10, \cr writePlots = TRUE, locusInfo = TRUE, path = "~/path/to/files")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
