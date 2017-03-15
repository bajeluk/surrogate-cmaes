#!/usr/bin/env Rscript
#
# A pairwise multi-comparison of algorithms with procedures described in
# Demsar (2006) and Garcia & Herrera (2008).
#
# Dependencies:
# * scmamp
# * Rgraphviz (for correct installation of the former)
# * optparse

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scmamp"))

# define command line arguments
option_list <- list(
                    make_option(c("-i", "--input"), default="algs.csv",
                                action="store",
                                help="Input file [default \"%default\"]"),
                    make_option(c("-o", "--output"), action="store",
                                help="Output file [default stdout]"),
                    make_option("--posthoc_test", action="store",
                                default="Friedman",
                                help=paste0("Posthoc test. One of [\"Friedman\", ",
                                            "\"Quade\", \"FriedmanAlignedRanks\"].",
                                            " [default \"Friedman\"]")),
                    make_option("--correction", action="store",
                                default="BergmannHommel",
                                help=paste0("P-values correction method. One of [",
                                            "\"BergmannHommel\", \"Shaffer\"].",
                                            " [default \"BergmannHommel\"]"))

               )

opt_parser = OptionParser(usage="usage: %prog [options]",
                          option_list=option_list, 
                          add_help_option=TRUE, prog=NULL, description="",
                          epilogue="")

# parse command line arguments
opts <- parse_args(opt_parser, print_help_and_exit=TRUE)

# validate arguments
if (!is.null(opts$output)) {
  out_fname <- opts$output
} else {
  # serves to write output to stdout
  out_fname <- ""
}

stopifnot(any(opts$posthoc_test == c("Friedman", "Quade", "FriedmanAlignedRanks")))
stopifnot(any(opts$correction == c("BergmannHommel", "Shaffer")))

if (opts$posthoc_test == "Friedman") {
  posthoc_test_fun = "friedmanPost"
} else if (opts$posthoc_test == "Quade") {
  posthoc_test_fun = "quadePost"
} else {
  posthoc_test_fun = "friedmanAlignedRanks"
}

adj_fun = paste0("adjust", opts$correction)

# read data
fvalues <- read.table(opts$input, header=TRUE, sep=",")

# perform the posthoc test
pv.matrix = do.call(posthoc_test_fun, list(fvalues))

# adjust p-values for multi-hypothesis testing
pv.adj = do.call(adj_fun, list(pv.matrix))

# output results
write.table(pv.matrix, out_fname, sep=",", quote=TRUE, eol="\n", row.names=TRUE,
            col.names=TRUE)
