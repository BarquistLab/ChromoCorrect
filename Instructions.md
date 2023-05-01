# Installation

1. Install R by going to the following website and following the installation instructions: https://cran.r-project.org/

2. Install RStudio by going to the following website and following the installation instructions: https://www.rstudio.com/products/rstudio/download/

3. Install the `devtools` package by running the following code in R:

install.packages("devtools")

4. Load the `devtools` package by running the following code:

```{r}
library(devtools)
```

5. Install the `ChromoCorrect` package from GitHub by running the following code:

```{r}
install_github("username/ChromoCorrect")
```

Note: Replace `username` with the actual username of the GitHub account that contains the `ChromoCorrect` package.

# Usage

1. Load the `ChromoCorrect` package by running the following code:

```{r}
library(ChromoCorrect)
```

2. To run the first script, navigate to the folder that contains your input files and run the following code:

```{r}
ChromoCorrect_script1("input_file1.csv", "input_file2.csv")
```

Note: Replace `input_file1.csv` and `input_file2.csv` with the actual names of your input files.

3. To run the second script, navigate to the folder that contains your input files and run the following code:

```{r}
ChromoCorrect_script2("input_file1.csv", "input_file2.csv")
```

Note: Replace `input_file1.csv` and `input_file2.csv` with the actual names of your input files.

4. To run the R Shiny app, run the following code:

```{r}
app()
```

This will launch the R Shiny app, which can be used to interactively analyze your data.

Note: Make sure that you have the required input files in the correct format in the same folder as the app file.

# Troubleshooting

If you encounter any errors while using the `ChromoCorrect` package, try reinstalling the package and restarting R. If the problem persists, please consult the package documentation or contact the package author for support.
