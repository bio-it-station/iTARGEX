# iTARGEX
This is git repo for iTARGEX (identification of Trait-Associated Regulatory Genes via mixture regression using EXpectation maximization). iTARGEX is a bioinformatics tool to help biological researchers for *in silico* prediction of novel regulators linked to biological traits. iTARGEX is mainly writen in R with embedded code in C for accelerating EM. In order to execute iTARGEX, the packages and its dependencies must be installed separately.


# Install
First, you need to install R enviroment in your computer. iTARGEX utlized Rcpp to accelerate the numerical calculation, so you need to install C++ compiler (e.g. gcc, clang). Next, you can install R packages and its dependencies using the following command.

```R
pkgs <- c("data.table", "docopt", "doParallel", "foreach", "mixtools", "Rcpp", "scales", "weights")
install.packages(pkgs)
```

# Quick Start
The simplest method to run iTARGEX is run `itargex.R` script without any options.

```sh
./bin/itargex.R trait.txt deletome.txt output_dir
```

This may take a few minutes to hours depending on your machine spec and data size. After execution, you can find the result under `output_dir`. The siginicant regulator can be found in `cor_sig_local.csv`. You can use any csv reader such as Excel or Libreoffice to open it. The first column is the name of signifiant regulator, followed by several columns.

| Column name | Description                                                    |
| ----------- | ---------------------------------------------------------------|
| cor         | Weighted pearson correlation coefficient                       |
| beta        | Slope of the significant component                             |
| pval_beta   | Significant level of beta (or cor)                             |
| adjp_beta   | Bonferroni-adjusted *p*-value                                  |
| ratio       | Proportion of datapoints estimated in the significant component|
| iter        | Iteration times during EM convergence                          |

# Advance usage
iTARGEX can be executed with serval options.

```
iTARGEX (Identify trail associated regulator via expectation–maximization).

Usage:
  itargex.R [options] <input-trait> <deletome-data> <output-folder>
  itargex.R em [--iter=<int> --ncpus=<int>] <input-trait> <deletome-data> <output-folder>
  itargex.R analysis [--ratio=<float> --ncpus=<int>] <input-trait> <deletome-data> <output-folder>
  itargex.R -h | --help
  itargex.R --version

Options:
  -h --help        Show this screen.
  --version        Show version.
  --ratio=<float>  Min ratio of regulator [default: 0.2].
  --iter=<int>     Max iteration time of EM algorithm [default: 1000].
  --ncpus=<int>    Number of CPU for parallel execution [default: 4].
```

If you want to execute iTARGEX for many times to search for the global optima of mixture models. You can run iTARGEX in two-step style to save the runtime.

```sh
# run EM step to esitmate the mixture weight of each datapoints
./bin/itargex.R em trait.txt deletome.txt output_dir
# run analyis step to identify significant regulators
./bin/itargex.R analysis --ratio 0.2 trait.txt deletome.txt output_dir
```

In addition, You can execute the second command for many time with different `--ratio`. Please notice that for any time the model runs it will overwrite the original result. We suggest you copy the ouput folder to other place to prevent overwrite.
