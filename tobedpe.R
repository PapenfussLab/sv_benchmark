# Exports the SV calls in a VCF file to BEDPE format
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
library(stringr)
library(VariantAnnotation)
library(devtools)
#library(optparse)
install_github("PapenfussLab/StructuralVariantAnnotation")
library(StructuralVariantAnnotation)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
	stop("Usage: Rscript tobedpe.R input.vcf output.bedpe")
}
vcf <- readVcf(args[1], "")
# Output BEDPE for use by circos
gr <- breakpointRanges(vcf)
bedpe <- data.frame(
    chrom1=seqnames(gr),
    start1=as.integer(start(gr) - 1),
    end1=as.integer(end(gr)),
    chrom2=seqnames(partner(gr)),
    start2=as.integer(start(partner(gr)) - 1),
    end2=as.integer(end(partner(gr))),
    name=names(gr),
    score=ifelse(is.nan(gr$QUAL) | is.na(gr$QUAL), 0, gr$QUAL),
    strand1=strand(gr),
    strand2=strand(partner(gr))
    )
write.table(bedpe, args[2], quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
