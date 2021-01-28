#!/usr/bin/env Rscript

"
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
  --ncpus=<int>    Number of CPU for parallel execution [default: 4].
" -> doc

library(docopt)
args <- docopt(doc, version = "0.1")

# casting to correct type
tryCatch(
    {
        args$iter <- as.integer(args$iter)
        args$ratio <- as.numeric(args$ratio)
        args$ncpus <- as.integer(args$ncpus)
    },
    warning = function(c) {
        message("Error: option should be numeric value.")
        quit(save = "no", status = 1)
    }
)

# validate options
if (args$ncpus < 1) {
    message("Error: --ncpus should larger than 1.")
    quit(save = "no", status = 1)
} else if (args$iter < 500) {
    message("Error: --iter should larger than 500.")
    quit(save = "no", status = 1)
} else if (args$ratio > 1 | args$ratio < 0) {
    message("Error: --ratio should between 0.0 ~ 1.0")
    quit(save = "no", status = 1)
}

orig_args <- commandArgs()
args0 <- normalizePath(dirname(strsplit(orig_args[sapply(orig_args, startsWith, prefix = "--file")], "=")[[1]][2]))

pos_args <- c(args$input_trait, args$deletome_data, args$output_folder)
if (args$em | !args$analysis) {
    cmd_args <- c(file.path(args0, "itargex_em.R"), "--iter", args$iter, "--ncpus", args$ncpus, pos_args)
    system2("Rscript", cmd_args)
}
if (args$analysis | !args$em) {
    cmd_args <- c(file.path(args0, "itargex_analysis.R"), "--ratio", args$ratio, "--ncpus", args$ncpus, pos_args)
    system2("Rscript", cmd_args)
}
