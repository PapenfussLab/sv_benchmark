source("sv_benchmark.R")
source("libplot.R")
library(rtracklayer)
library(dplyr)
library(openxlsx)
library(stringr)

maxgap <- 180
ignore.strand <- FALSE
sizemargin <- 0.25
maxSpanningFragmentSize <- 500
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
	gr$insLen <- 0
	gr$svLen <- NA_integer_
	gr$vcfId <- paste0(sheetname, "_", seq_along(gr))
	return(gr)
})
names(publishedgrs) <- str_replace(drsheetnames, stringr::fixed(" (DR)"), "")
publishedgrs <- sapply(names(publishedgrs), function(sample) {
	gr <- publishedgrs[[sample]]
	seqlevelsStyle(gr) <- "UCSC"
	#gr <- subsetbed(gr, publishedcgrs[[sample]], maxgap)
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)

#########################
# Load VCFs
metadata <- LoadMetadata(paste0(rootdir, "data.neo"))
metadata$samplename <- metadata$CX_SAMPLE %na% toupper(str_match(metadata$CX_BAM, "([^/]+).s..bam$")[,2] %na%
  # hack to match on fragment size
  (metadata %>% left_join(metadata %>% filter(Id %in% c("778", "got3", "t1000")), by="CX_READ_FRAGMENT_LENGTH"))$Id.y)
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.neo"), metadata=metadata, existingList=vcfs)
vcfs <- sapply(names(vcfs), function(id) {
	sample <- (metadata %>% filter(Id==id))$samplename
	gr <- vcfs[[id]]
	seqlevelsStyle(gr) <- "UCSC"
	#gr <- subsetbed(gr, publishedcgrs[[sample]], maxgap)
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)#########################
# Match calls
.summarymatches <- function(id, includeFiltered) {
	write(paste0("Processing ", id), stderr())
	sample <- (metadata %>% filter(Id==id))$samplename
	gr <- vcfs[[id]]
	gr$incgr <- overlapsAny(gr, publishedcgrs[[sample]], maxgap=1000)
	gr$bothgr <- gr$incgr + partner(gr)$incgr
	pubgr <- publishedgrs[[sample]]
	gridssgr <- vcfs[[(metadata %>% filter(samplename==sample & CX_ALIGNER == "bwamem" & StripCallerVersion(CX_CALLER) == "gridss"))$Id]]
	gridssgr <- gridssgr[gridssgr$FILTER %in% c(".", "PASS")]
	calls <- ScoreVariantsFromTruthVCF(gr, pubgr, includeFiltered=includeFiltered, maxgap, ignore.strand, sizemargin)
	callsgridss <- ScoreVariantsFromTruthVCF(gr, gridssgr, includeFiltered=includeFiltered, maxgap, ignore.strand, sizemargin)

	#spanninghits <- findSpanningBreakpoints(gr, pubgr, maxgap, ignore.strand, sizemargin, maxSpanningFragmentSize=maxSpanningFragmentSize, matchDirection=FALSE)
	#calls$truth$spanningtp <- FALSE
	#calls$truth$spanningtp[spanninghits$squeryHits] <- TRUE
	#calls$truth$tp <- calls$truth$tp | calls$truth$spanningtp
	
	# filter small events that did not match a published call
	# filter false positives outside of CGRs
	misses <- !calls$calls$tp & !callsgridss$calls$tp & gr$bothgr & (is.na(gr$svLen) | abs(gr$svLen) >= minimumEventSize)
	
	callsummary <- data.frame(
	  Id=id,
	  pubmatches=sum(calls$truth$tp) / 2,
	  gridssmatches=sum(callsgridss$truth$tp) / 2,
	  misses=sum(misses) / 2,
	  pubcount=length(pubgr) / 2,
	  gridsscount=length(gridssgr) / 2
	  #spanning=nrow(spanninghits) / 2,
	  #spanningMeanFragmentSize=mean(spanninghits$fragmentSize)
	)
	return(callsummary)
}
summarylist <- lapply((metadata %>% filter(!is.na(CX_CALLER) & Id %in% names(vcfs)))$Id, function(id) .summarymatches(id, TRUE))
summarydf <- rbind_all(summarylist) %>%
  left_join(metadata %>% select(Id, CX_CALLER, CX_ALIGNER)) %>%
  mutate(caller=StripCallerVersion(CX_CALLER))
summary <- summarydf %>% group_by(CX_CALLER, CX_ALIGNER) %>%
  summarise(pubsens=sum(pubmatches)/sum(pubcount), gridsssens=sum(gridssmatches)/sum(gridsscount), misses=sum(misses))

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


