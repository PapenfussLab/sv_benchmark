source("sv_benchmark.R")
source("libplot.R")
library(GenomicRanges)
library(dplyr)
library(stringr)
library(xlsx)

rnaxlsx <- read.xlsx(paste0(rootdir, "input.HCC1395/INTEGRATE/Supplemental_Table2.xlsx"), 1, startRow=6)
rnagr5p <- GRanges(
	seqnames=rnaxlsx$X5p_chr,
	ranges=IRanges(start=rnaxlsx$X5p_fusion_junction, width=1),
	strand=rnaxlsx$X5p_strand)
rnagr3p <- GRanges(
	seqnames=rnaxlsx$X3p_chr,
	ranges=IRanges(start=rnaxlsx$X3p_fusion_junction, width=1),
	strand=rnaxlsx$X3p_strand)
rnagr <- c(rnagr5p, rnagr3p)
rnagr$partner <- c(paste0("three", seq_along(rnagr5p)), paste0("five", seq_along(rnagr5p)))
names(rnagr) <- c(paste0("five", seq_along(rnagr5p)), paste0("three", seq_along(rnagr5p)))
rnagr$svLen <- NA
rnagr$insLen <- NA
rnagr$vcfId <- NA
seqlevelsStyle(rnagr) <- "UCSC"

dnaxlsx <- read.xlsx(paste0(rootdir, "input.HCC1395/INTEGRATE/Supplemental_Table3.xlsx"), 1, startRow=6)
dnagr5p <- GRanges(
	seqnames=dnaxlsx$chr_5p,
	ranges=IRanges(start=dnaxlsx$bk_5p, width=1))
dnagr3p <- GRanges(
	seqnames=dnaxlsx$chr_3p,
	ranges=IRanges(start=dnaxlsx$bk_3p, width=1))
dnagr <- c(dnagr5p, dnagr3p)
dnagr$partner <- c(paste0("three", seq_along(dnagr5p)), paste0("five", seq_along(dnagr5p)))
names(dnagr) <- c(paste0("five", seq_along(dnagr5p)), paste0("three", seq_along(dnagr5p)))
dnagr$svLen <- NA
dnagr$insLen <- NA
dnagr$vcfId <- NA
seqlevelsStyle(dnagr) <- "UCSC"


dnagr$svLen <- abs(start(dnagr) - start(partner(dnagr)))
rnagr$svLen <- abs(start(rnagr) - start(partner(rnagr)))

#dnagr <- dnagr[dnagr$svLen >= 1000000]
#rnagr <- rnagr[dnagr$svLen >= 1000000]
#Problem 1: RNASeq calls are at exon boundaries
#DNASeq calls are at genomic breakpoints - these can be way off
# eg: UNC5C-STPG2 positions differ by over 202kpb


matches <- findBreakpointOverlaps(rnagr, dnagr, maxgap=205000)
rnagr$dnaindex <- NA
rnagr$dnaindex[matches$queryHits] <- matches$queryHits
dnagr$rnaindex <- NA
dnagr$rnaindex[matches$subjectHits] <- matches$subjectHits


metadata <- LoadMetadata(paste0(rootdir, "data.HCC1395"))
metadata <- metadata %>% filter(str_detect(CX_BAM, "bwa.sc.bam"))
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.HCC1395"), metadata=metadata)
vcf <- vcfs[[1]]
vcf <- vcf[is.na(vcf$svLen) | abs(vcf$svLen) > 100,]
# convert coordinates

dnacalls <- ScoreVariantsFromTruthVCF(vcf, dnagr, TRUE, maxgap=1000, ignore.strand=TRUE)
rnacalls <- ScoreVariantsFromTruthVCF(vcf, rnagr, TRUE, maxgap=250000, ignore.strand=TRUE)
dnagr$gridss <- NA
dnagr$gridss <- dnacalls$truth$QUAL
rnagr$gridss <- NA
rnagr$gridss <- rnacalls$truth$QUAL

# Found by GRIDSS and INTEGRATE
table(rnagr$gridss > 0, !is.na(rnagr$dnaindex))
