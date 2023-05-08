# Installation
1. Install R and RStudio by going to this website and following the installation instructions: https://www.rstudio.com/products/rstudio/download/
2. Install and load the `devtools` package by running the following code in R:
```{r}
install.packages("devtools")
library(devtools)
```
3. Install the `ChromoCorrect` package from GitHub by running the following code:
```{r}
install_github("gerisullivan/ChromoCorrect")
```

## Folder organisation
To organise your input and output files, we recommend creating a folder on your local computer. Within this folder, create two sub-folders: **logFCs**, and **readcounts**.

In **logFCs**, put your output files either from Bio::TraDIS or other analysis pipelines. This requires the columns *locus_tag* and *logFC*, with loci in the same order as they are on the chromosome (if you are using default TraDIS output files, the order is correct).

In **readcounts**, put your traDIS insert site files, or files containing read counts from other analysis pipelines that require correction. The files require a *locus_tag* column and a *read_count* column for the script to run. The script also assumes the loci are in the same order as they are on the chromosome.  

## Launching the app
1. Load the `ChromoCorrect` package by running the following code:
```{r}
library(ChromoCorrect)
```
2. Set your working directory so you can locate your output:
```{r}
setwd("/path/to/ChromoCorrectdir")
```

3. Launch the app with the aptly named command:
```{r}
launch_app()
```
This will launch the R Shiny app, which can be used to interactively analyse your data.

## Running the normalisation locally
1. Load the `ChromoCorrect` package by running the following code:
```{r}
library(ChromoCorrect)
```
2. Set your working directory so you can locate your output:
```{r}
setwd("/path/to/ChromoCorrectdir")
```

3. Diagnose your output files for chromosomal location bias:
```{r}
# assuming your logFC files are in a folder called logFCs, the default will work:
detect_bias()
```
A file called 'detect_bias.png' will be saved in your logFCs folder. Look at it to determine which, if any, of your conditions are affected by chromosomal location bias.

4. Correct your chromosomal location bias, part 1:

This assumes you have at least two conditions of read counts per run of the script, with replicate information in the form of _1, _2 before the suffix.
e.g. MH_1.tradis_gene_insert_sites.csv, MH_2.tradis_gene_insert_sites.csv,
Cip_1.tradis_gene_insert_sites.csv, Cip_2.tradis_gene_insert_sites.csv.
```{r}
# import your read counts into one file:
readcounts <- structure_rc(csvpath = "/readcounts", getLocusInfo = TRUE, suffix = ".tradis_gene_insert_sites.csv")
```
This will save the locus information as 'locusInfo.tsv' for the next script.

5. Correct your chromosomal location bis, part 2:

```{r}
# normalise your read counts:
norm_readcounts <- normalise_bias(readcounts)
# It will prompt you to choose which of your read count files is your control.
```

# Example data
The ciprofloxacin data set published with this paper is available under /ChromoCorrect/inputData/. It can be used within the app and locally.

# Troubleshooting
If you encounter any errors while using the `ChromoCorrect` package, try reinstalling the package and restarting R. If the problem persists, please consult the package documentation or contact the package author for support.
