source("config.R")
source("sv_benchmark.R")
source("libplot.R")
source("shinyCache.R")
library(shiny)
library(stringr)

enableBookmarking(store = "url")
options("R.cache::compress" = TRUE)

# load all metadata
md <- lapply(list.files(dataLocation, pattern = "data.*", full.names = TRUE), function(dir) LoadCachedMetadata(dir))
names(md) <- str_replace(list.files(dataLocation, pattern = "data.*"), "data.", "")
# hack fix for:
# Error: incompatible type (data index: 4, column: 'CX_MULTIMAPPING_LOCATIONS', was collecting: integer (dplyr::Collecter_Impl<13>), incompatible with data of type: character
md <- lapply(md, function(x) { x$CX_MULTIMAPPING_LOCATIONS <- as.integer(x$CX_MULTIMAPPING_LOCATIONS) ; x})

# load blacklist bed
if (!exists("lrblacklistgr")) {
  write("Loading blacklist", stderr())
  lrblacklistgr <- list(
    None=GRanges(),
    DAC=import(paste0(dataLocation, "/input.common/wgEncodeHg19ConsensusSignalArtifactRegions.bed")),
    Duke=import(paste0(dataLocation, "/input.common/wgEncodeDukeMapabilityRegionsExcludable.bed"))
  )
	seqlevelsStyle(lrblacklistgr$DAC) <- "UCSC"
	seqlevelsStyle(lrblacklistgr$Duke) <- "UCSC"
}

# load repeatmasker annotations
if (!exists("grrm")) {
	key = list(repeatmaskermergedfile=paste0(dataLocation, "/input.common/repeatmasker-hg19.fa.out.gz"))
	grrm <- loadCache(key=key)
	if (is.null(grrm)) {
	  write(paste("Loading ", key$repeatmaskermergedfile), stderr())
		grrm <- import.repeatmasker.fa.out(key$repeatmaskermergedfile)
		saveCache(grrm, key=key)
	}
}

# helper functions
PrettyAligner <- function(dataaligners) {
    dataaligners <- unique(dataaligners)
    dataaligners <- dataaligners[!is.na(dataaligners)]
    ka <- knownaligners[knownaligners %in% dataaligners]
    return(c("Most sensitive" = "best", ka, "N/A" = ""))
}
withnames <- function(v, n) { names(v) <- n; return(v) }
.primaryHumanOnly <- function(gr, metadata) {
	if (length(gr) > 0) {
		seqlevelsStyle(gr) <- "UCSC"
		# filter to primary chromosomes
		gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")) & seqnames(partner(gr)) %in% paste0("chr", c(1:22, "X", "Y")),]
	}
	repeatHits <- findOverlaps(gr, grrm, select="first")
	gr$repeatClass <- ifelse(is.na(repeatHits), "", grrm[repeatHits %na% 1]$repeatClass)
    return(gr)
}
.primaryHumanOnly_blacklist <- function(gr, metadata, blacklist) {
  blacklistgr <- lrblacklistgr[[blacklist]]
	gr <- .primaryHumanOnly(gr, metadata)
	gr <- gr[!overlapsAny(gr, blacklistgr) & !overlapsAny(partner(gr), blacklistgr),]
  return(gr)
}

# callers that have data points for all simulation conditions and NA12878
# TODO: replace this with a computed field that tells us what
# data set each caller is missing
fulldatacallers <- c(
	"breakdancer",
	"cortex",
	"crest",
	"delly",
	"gasv",
	"gridss",
	"hydra",
	"lumpy",
	"manta",
	"pindel",
	"socrates"
)
eventtypes = c("Deletion"="DEL", "Insertion"="INS", "Inversion"="INV", "Tandem Duplication"="DUP") # "BP"
knownaligners <- c("bowtie2", "bwa mem"="bwamem", "novoalign")
simfacets <- c(
				"Read Length"="CX_READ_LENGTH",
				"Read Depth"="CX_READ_DEPTH",
				"Fragment Size"="CX_READ_FRAGMENT_LENGTH",
				"Aligner"="aligner",
                "Call Set" = "CallSet")
dataoptions <- list()
dataoptions$maxgap <- 200
dataoptions$ignore.strand <- TRUE
dataoptions$sizemargin <- 0.25
dataoptions$ignore.duplicates <- FALSE
dataoptions$ignore.interchromosomal <- TRUE
dataoptions$mineventsize <-51
dataoptions$maxeventsize <- NULL
dataoptions$requiredHits <- 1
dataoptions$datadir <- names(md)
dataoptions$grtransform <- .primaryHumanOnly
simoptions <- dataoptions
simoptions$datadir <- c("Read Depth" = "rd", "Read Length" = "rl", "Fragment Size" = "fs")
simoptions$datasetslicecol <- c("rd" = "CX_READ_DEPTH", "rl" = "CX_READ_LENGTH", "fs" = "CX_READ_FRAGMENT_LENGTH")
simoptions$mineventsize <- c(0, 51)
# simoptions$grtransform <- # filter small events calls on breakpoint data sets to remove spurious indels caused by sequence homology around breakpoints?
lroptions <- dataoptions
lroptions$requiredHits <- 3
lroptions$datadir <- lroptions$datadir[!(lroptions$datadir %in% simoptions$datadir)]
# lapply doesn't quite work since it doesn't play nicely with R.cache
#lroptions$grtransform <- lapply(names(lrblacklistgr), function(blacklist) function(gr, metadata) .primaryHumanOnly_blacklist(gr, metadata, blacklist))
lroptions$grtransform <- list(
  None=function(gr, metadata) .primaryHumanOnly_blacklist(gr, metadata, "None"),
  DAC=function(gr, metadata) .primaryHumanOnly_blacklist(gr, metadata, "DAC"),
  Duke=function(gr, metadata) .primaryHumanOnly_blacklist(gr, metadata, "Duke"))
                                
                                