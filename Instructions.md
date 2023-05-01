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
# Launching the app
1. Load the `ChromoCorrect` package by running the following code:
```{r}
library(ChromoCorrect)
```
2. Launch the app with the aptly named command
```{r}
launch_app()
```
This will launch the R Shiny app, which can be used to interactively analyse your data.

# Example data
The ciprofloxacin data set published with this paper is available under /ChromoCorrect/data/. It can be used within the app and locally.

# Running the normalisation locally
## Folder organisation
To organise your input and output files, we recommend creating a folder and R Project specifically for the ChromoCorrect package. Within this folder, create two sub-folders: ‘logFCs’, and ‘readcounts’.
In ‘logFCs’, put your output files either from Bio::TraDIS or other analysis pipelines. This requires the columns *locus_tag* and *logFC*, with loci in the same order as they are on the chromosome. If you are using default TraDIS output files, the order is correct.
In ‘readcounts’, put your tradis insert site files or files containing read counts from other analysis pipelines. The files require a *locus_tag* column and a *read_count* column for the script to run. The script also assumes the loci are in the same order as they are on the chromosome.
## Run the normalisation
1. To make a new R Project, open an RStudio session and navigate to File --> New Project… --> Existing Directory and browse for your ChromoCorrect folder. Tick the box to open in a new session.

Place your input files in the "input" folder and your output files in the "output" folder. Any additional data files needed for the package should be placed in the "data" folder. This organization will make it easier for users to keep track of their input and output files, as well as any data files that are needed for the package.
1. Load the `ChromoCorrect` package by running the following code:
```{r}
library(ChromoCorrect)
```
2. If you have Bio::TraDIS output files
```{r}
ChromoCorrect_script1("input_file1.csv", "input_file2.csv")
```

3. To run the second script, navigate to the folder that contains your input files and run the following code:

```{r}
ChromoCorrect_script2("input_file1.csv", "input_file2.csv")
```
# Troubleshooting
If you encounter any errors while using the `ChromoCorrect` package, try reinstalling the package and restarting R. If the problem persists, please consult the package documentation or contact the package author for support.
