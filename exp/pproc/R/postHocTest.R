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
                    make_option(c("-p", "--padj_out"), action="store",
                                default="",
                                help=paste0("Output file for adjusted p-vals ",
                                            "(default stdout).")),
                    make_option(c("-s", "--summary_out"), action="store",
                                default="",
                                help=paste0("Output file for the test summary ",
                                            "(default stdout).")),
                    make_option(c("-t", "--posthoc_test"), action="store",
                                default="friedman",
                                help=paste0("Posthoc test. One of ",
                                            "\"friedman\" (default), ",
                                            "\"wilcoxon\", ",
                                            "\"quade\", ",
                                            "\"aligned ranks\".")),
                    make_option(c("-c", "--correction"), action="store",
                                default="bergmann",
                                help=paste0("P-values correction method. One of ",
                                            "\"bergmann\" (default), ",
                                            "\"shaffer\", or \"bonferrroni\"."))

               )

opt_parser = OptionParser(usage="usage: %prog [options]",
                          option_list=option_list, 
                          add_help_option=TRUE, prog=NULL, description="",
                          epilogue="")

# parse command line arguments
opts <- parse_args(opt_parser, print_help_and_exit=TRUE)

stopifnot(any(opts$posthoc_test == c("friedman", "quade",
                                     "wilcoxon", "aligned ranks")))
stopifnot(any(opts$correction == c("bergmann", "shaffer", "bonferroni")))

# read data
fvalues <- read.table(opts$input, header=FALSE, sep=",")
fvalues <- (max(fvalues)-fvalues)/(max(fvalues)-min(fvalues))

# perform the posthoc test
res <- postHocTest(data=fvalues, test=opts$posthoc_test,
                   control=NULL, use.rank=TRUE,
                   correct=opts$correction)

# output adjusted p-values
write.table(res$corrected.pval, opts$padj_out, na="NaN", sep=",",
            quote=FALSE, eol="\n", row.names=FALSE, col.names=FALSE)

# output summary
write.table(res$summary, opts$summary_out, na="NaN", sep=",",
            quote=FALSE, eol="\n", row.names=FALSE, col.names=FALSE)
