library(devtools)
library(testthat)
library(roxygen2)
library(assertthat)
install_github("d-cameron/StructuralVariantAnnotation")
library(StructuralVariantAnnotation)
library(GenomicRanges)
library(reshape)
library(stringr)
library(dplyr)

minsize <- 1000
vcf <- readVcf("C:/dev/somatic/MELA0248.gridss.vcf", "hg19")
gr <- breakpointRanges(vcf)
# strip out small events
vcf <- vcf[is.na(gr$svLen) | abs(gr$svLen > minsize),]
vcf <- vcf[rowRanges(vcf)$FILTER %in% c(".", "PASS"),]

df <- unpack(info(vcf))
vcf <- vcf[df$SR.1 + df$RSR.1 + df$RP.1 + df$ASSR.1 + df$ASRP.1 == 0,]

# Convert to BEDPE
gr <- breakpointRanges(vcf)
bedpe <- data.frame(
	chrom1=seqnames(gr),
	start1=start(gr),
	end1=end(gr),
	chrom2=seqnames(partner(gr)),
	start2=start(partner(gr)),
	end2=end(partner(gr)),
	name=names(gr),
	score=gr$QUAL,
	strand1=strand(gr),
	strand2=strand(partner(gr))
	)
bedpe <- bedpe[str_detect(bedpe$name, "gridss[0-9]+o"),] # Just the lower of the two breakends
write.table(bedpe, "C:/dev/somatic/MELA0248.gridss.hq.somatic.bedpe", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
