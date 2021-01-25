# iTARGEX
This is git repo for iTARGEX. iTARGEX is a program help biological researcher to identify trail associated regulator. iTARGEX is mainly writen in R and depend on some packages. In order to execute iTARGEX you need install R and the its dependency.

# Install
First, you need to install R enviroment in your computer. iTARGEX utlized Rcpp to accelerate the numerical calculation, so you need to install C++ compiler (e.g. gcc, clang). Next, you can install R package dependency with the following command.

```R
pkgs <- c("data.table", "docopt", "doParallel", "foreach", "mixtools", "Rcpp", "scales")
install.packages(pkgs)
```

# Quick Start
The simplest method to run iTARGEX is run `itargex.R` script without any options.

```sh
./bin/itargex.R trait.txt deletome.txt output_dir
```

This may take a few minites to hours depend on your machine spec and data size. After execution, you can find the result under `output_dir`. The siginicant regulator can be found in `cor_sig_local.csv`. You can use any csv reader like Excel or Libreoffice to open it. The first column is the name of signifiant regulator, followed by several column.

| Column name | Description                              |
| ----------- | ---------------------------------------- |
| cor         | weighted pearson correlation coefficient |
| beta        | slope of significant component           |
| pval_beta   | significant level of beta (or cor)       |
| adjp_beta   | Bonferroni corrected *p*-value           |
| ratio       | ratio of regulated gene                  |
| iter        | iteration number during EM algorithm     |

# Advance usage
iTARGEX can execute with serval options.

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
  --iter=<int>     Max iteration number of EM algorithm [default: 1000].
  --ncpus=<int>    Number of CPU for parallel execuation [default: 4].
```

If you wish to execute iTARGEX for many times to fine tune the min ratio cutoff (`--ratio`). You can run iTARGEX in two step style to save the runtime.

```sh
# run EM step to esitmate the weight
./bin/itargex.R em trait.txt deletome.txt output_dir
# run analyis step to identify significant regulator
./bin/itargex.R analysis --ratio 0.2 trait.txt deletome.txt output_dir
```

You can execute the second command for many time with different `--ratio`. Please notice that run multiple iTARGEX will overwrite the orignial result. We suggest you copy the ouput folder to other place to prevent overwrite.