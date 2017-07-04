#!/usr/bin/env Rscript
library("optparse")
parser <- OptionParser(option_list=list(
	make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate plots")
))
arguments = parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args
write(commandArgs(trailingOnly=TRUE), stderr())
write(names(opt), stderr())
write(opt$plot, stderr())
