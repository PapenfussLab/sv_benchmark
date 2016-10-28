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

# helper functions
PrettyAligner <- function(simdatadirs) {
    dataaligners <- unique(md[[ds]]$CX_ALIGNER)
    dataaligners <- dataaligners[!is.na(dataaligners)]
    ka <- knownaligners[knownaligners %in% dataaligners]
    return(c("Most sensitive" = "best", ka, "N/A" = ""))
}
withnames <- function(v, n) { names(v) <- n; return(v) }

# callers that have data points for all simulation conditions and NA12878
# TODO: replace this with a computed field that tells us what
# data set each caller is missing
fulldatacallers <- c(
	"breakdancer",
	"cortex",
	"crest",
	"delly",
	"gridss",
	"hydra",
	"lumpy",
	"manta",
	"pindel",
	"socrates"
)
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
dataoptions$ignore.duplicates <- TRUE
dataoptions$ignore.interchromosomal <- TRUE
dataoptions$mineventsize <-51
dataoptions$maxeventsize <- NULL
dataoptions$requiredHits <- 1
dataoptions$datadir <- names(md)
dataoptions$vcftransform <- function(id, metadata, vcfs) {
    gr <- vcfs[[id]]
    if (!is.na((metadata %>% filter(Id == id))$CX_CALLER)) {
        # only looking at intrachromosomal calls
        gr <- gr[!is.na(gr$svLen),]
    }
    return(gr)
}
simoptions <- dataoptions
simoptions$datadir <- c("Read Depth" = "rd", "Read Length" = "rl", "Fragment Size" = "fs")
simoptions$datasetslicecol <- c("rd" = "CX_READ_DEPTH", "rl" = "CX_READ_LENGTH", "fs" = "CX_READ_FRAGMENT_LENGTH")
simoptions$mineventsize <- c(0, 51)
simoptions$vcftransform <- function(id, metadata, vcfs) {
    gr <- dataoptions$vcftransform(id, metadata, vcfs)
    if (!is.na((metadata %>% filter(Id == id))$CX_CALLER)) {
        #if (str_detect((metadata %>% filter(Id == id))$CX_REFERENCE_VCF_VARIANTS, "BP")) {
            # filtering small events calls to remove spurious indels caused by sequence homology around breakpoints
        #	gr <- gr[abs(gr$svLen) >= 50,]
        #}
    }
    return(gr)
}
lroptions <- dataoptions
lroptions$requiredHits <- 3
lroptions$datadir <- lroptions$datadir[!(lroptions$datadir %in% simoptions$datadir)]

