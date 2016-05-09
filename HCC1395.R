source("sv_benchmark.R")
source("libplot.R")
library(GenomicRanges)
library(dplyr)
library(stringr)
library(xlsx)

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/projects/sv_benchmark/", "~/projects/sv_benchmark/")

fusionxlsx <- read.xlsx(paste0(rootdir, "input.HCC1395/INTEGRATE.Supplemental_Table2.xlsx"), 1)
fusiongr5p <- GRanges(
	seqnames=fusionxlsx$X5p_chr,
	ranges=IRanges(start=fusionxlsx$X5p_fusion_junction, width=1),
	strand=fusionxlsx$X5p_strand)
fusiongr3p <- GRanges(
	seqnames=fusionxlsx$X3p_chr,
	ranges=IRanges(start=fusionxlsx$X3p_fusion_junction, width=1),
	strand=fusionxlsx$X3p_strand)
fusiongr <- c(fusiongr5p, fusiongr3p)
fusiongr$partner <- c(paste0("three", seq_along(fusiongr5p)), paste0("five", seq_along(fusiongr5p)))
names(fusiongr) <- c(paste0("five", seq_along(fusiongr5p)), paste0("three", seq_along(fusiongr5p)))
