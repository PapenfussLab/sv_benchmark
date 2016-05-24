source("sv_benchmark.R")
source("libplot.R")
library(rtracklayer)
library(dplyr)
library(openxlsx)
library(stringr)

maxgap <- 180
ignore.strand <- FALSE
sizemargin <- 0.25
maxSpannnedFragmentSize <- 375
minimumEventSize <- 500

neochromosome_supp3 <- paste0(rootdir, "input.neo/mmc3.xlsx")
neochromosome_supp4 <- paste0(rootdir, "input.neo/mmc4.xlsx")
#########################
# Extract published CGRs
cgrsheetnames <- getSheetNames(neochromosome_supp4)
cgrsheetnames <- cgrsheetnames[str_detect(cgrsheetnames, stringr::fixed("_CGRs"))]
publishedcgrs <- lapply(cgrsheetnames, function(sheetname) {
	dt <- read.xlsx(neochromosome_supp4, sheetname)
	gr <- GRanges(seqnames=dt$Chromosome, ranges=IRanges(start=dt$Start, end=dt$End))
	names(gr) <- dt$CGR.name
	return(gr)
})
names(publishedcgrs) <- str_replace(cgrsheetnames, stringr::fixed("_CGRs"), "")
#########################
# Extract published CNs
cnsheetnames <- getSheetNames(neochromosome_supp4)
cnsheetnames <- cnsheetnames[str_detect(cnsheetnames, stringr::fixed("_CN"))]
publishedcns <- lapply(cnsheetnames, function(sheetname) {
	dt <- read.xlsx(neochromosome_supp4, sheetname)
	gr <- GRanges(seqnames=dt$Chromosome,
		ranges=IRanges(start=dt$Start, end=dt$End),
		cgrname=dt$CGR.name,
		cn=dt$Copy.number)
	return(gr)
})
names(publishedcns) <- str_replace(cnsheetnames, stringr::fixed("_CN"), "")
#########################
# Extract published calls
drsheetnames <- getSheetNames(neochromosome_supp3)
drsheetnames <- drsheetnames[str_detect(drsheetnames, stringr::fixed("(DR)"))]
publishedgrs <- lapply(drsheetnames, function(sheetname) {
	dt <- read.xlsx(neochromosome_supp3, sheetname)
	gr <- GRanges(seqnames=c(dt$chrom1, dt$chrom2),
                ranges=IRanges(start=c(dt$start1, dt$start2), width=1),
                strand=c(as.character(dt$strand1), as.character(dt$strand2)),
                partner=c(paste0("row", seq_along(dt$chrom1), "_bp2"), paste0("row", seq_along(dt$chrom2), "_bp1")),
                QUAL=dt$nreads)
	names(gr) <- c(paste0("row", seq_along(dt$chrom1), "_bp1"), paste0("row", seq_along(dt$chrom2), "_bp2"))
	return(gr)
})
names(publishedgrs) <- str_replace(drsheetnames, stringr::fixed(" (DR)"), "")
publishedgrs <- sapply(names(publishedgrs), function(sample) {
	gr <- publishedgrs[[sample]]
	seqlevelsStyle(gr) <- "UCSC"
	gr <- subsetbed(gr, publishedcgrs[sample])
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)

#########################
# Load VCFs
metadata <- LoadMetadata(paste0(rootdir, "data.neo"))
metadata$samplename <- toupper(str_extract(md$CX_REFERENCE_BAM_TODO), "[^/]+.sc.bam")
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.neo"), metadata=metadata)
vcfs <- sapply(names(vcfs), function(id) {
	sample <- (metadata %>% filter(Id==id))$samplename
	gr <- vcfs[[id]]
	seqlevelsStyle(gr) <- "UCSC"
	gr <- subsetbed(gr, publishedcgrs[sample])
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)

#########################
# Match calls
lapply((metadata %>% filter(samplename==sample))$Id, function(id) {
	sample <- (metadata %>% filter(Id==id))$samplename
	gr <- vcfs[[id]]
	pubgr <- publishedgrs[sample]
	gridssgr <- vcfs[[(metadata %>% filter(samplename==sample & StripCallerVersion(CX_CALLER) == "gridss"))$Id]]
	hitpub <- ScoreVariantsFromTruthVCF(gr, pubgr, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin)
	hitgridss <- ScoreVariantsFromTruth(vcfs, gridssgr, includeFiltered=TRUE, maxgap, ignore.strand, sizemargin, truthgr=pubgr)

	calls <- hitpub$calls
	calls$publishedtp <- calls$tp
	calls$gridsstp <- hitgridss$calls$tp

	# TODO: look for pub spanning calls
	# 1) find overlapping breakends
	# 2) find call pairs such that query breakends Q1, Q2 have matches Q1=A1, Q2=B2 and A2<->B1 < maxSpannnedFragmentSize

	# TODO: filter small events that did not match pub set

	return (data.frame(
		Id=id,
		pubmatches=sum(calls$publishedtp),
		gridssmatches=sum(!calls$publishedtp & calls$gridsstp),
		misses=sum(!calls$gridsstp & !calls$publishedtp)
	))
})


# Neochromosome venn diagram

# Venn Diagram for each other caller

# list of overlaps

# Published, GRIDSS, ...
# 1
# 2
# 3
# 4
# 5
# 6
# 7
# 8
# High Confidence calls

# within main text: compound breakpoints identified by caller
# Published, GRIDSS, ...
# identify calls that could be transitive
# identify calls that are transitive


