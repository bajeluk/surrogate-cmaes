#!/usr/bin/env Rscript
#
# A pairwise multi-comparison of algorithms with procedures described in
# Demsar (2006) and Garcia & Herrera (2008).
#
# Dependencies:
# * scmamp
# * Rgraphviz (for correct installation of the former)
# * optparse
#
# How to install dependencies in an R session:
# > install.packages("optparse")
# > # moved from CRAN to bioconductor
# > source("http://www.bioconductor.org/biocLite.R")
# > biocLite("Rgraphviz")
# > install.packages("scmamp")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scmamp"))

# define command line arguments
option_list <- list(
                    make_option(c("-i", "--input"), default="algs.csv",
                                action="store",
                                help="Input file (default \"%default\")"),
                    make_option(c("-p", "--p_out"), action="store",
                                default="",
                                help=paste0("Output for the p-value ",
                                            "(default stdout).")),
                    make_option(c("-s", "--statistic_out"), action="store",
                                default="",
                                help=paste0("Output file for the statistic ",
                                            "(default stdout).")),
                    make_option(c("-t", "--test"), action="store",
                                default="friedman",
                                help=paste0("Multiple comparison test. Either ",
                                            "\"iman\" (default) ",
                                            "or \"friedman\"."))
               )

opt_parser = OptionParser(usage="usage: %prog [options]",
                          option_list=option_list, 
                          add_help_option=TRUE, prog=NULL, description="",
                          epilogue="")

# parse command line arguments
opts <- parse_args(opt_parser, print_help_and_exit=TRUE)

stopifnot(any(opts$test == c("friedman", "iman")))

# read data
fvalues <- read.table(opts$input, header=FALSE, sep=",")
fvalues <- (max(fvalues)-fvalues)/(max(fvalues)-min(fvalues))

# perform the posthoc test
res <- multipleComparisonTest(data=fvalues, test=opts$test)

# output the p-value
write(res$p.value, opts$p_out)

# output the statistic
write(res$statistic, opts$statistic_out)
