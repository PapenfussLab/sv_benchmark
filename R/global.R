source("sv_benchmark.R")
source("shinyCache.R")

knownaligners <- c("bowtie2", "bwa mem"="bwamem", "novoalign")
ds <- c("Read Depth"="rd", "Read Length"="rl", "Fragment Size"="fs")
md <- lapply(ds, function (data) LoadCachedMetadata(paste0("../data.", data)))
names(md) <- ds
withnames <- function(v, n) { names(v) <- n; return(v) }
PrettyAligner <- function(ds) {
	dataaligners <- unique(md[[ds]]$CX_ALIGNER)
	dataaligners <- dataaligners[!is.na(dataaligners)]
	ka <- knownaligners[knownaligners %in% dataaligners]
	return (c("Most sensitive"="best", ka, "N/A"=""))
}
facets <- c(
				"Read Length"="CX_READ_LENGTH",
				"Read Depth"="CX_READ_DEPTH",
				"Fragment Size"="CX_READ_FRAGMENT_LENGTH",
				"Aligner"="aligner",
				"Call Set"="CallSet")

simoptions <- list()
simoptions$ds <- ds
simoptions$maxgap <- 200
simoptions$ignore.strand <- TRUE
simoptions$sizemargin <- 0.25
simoptions$ignore.duplicates <- TRUE
simoptions$ignore.interchromosomal <- TRUE
simoptions$mineventsize <- c(0, 51)
simoptions$maxeventsize <- NULL
simoptions$vcftransform <- function(id, metadata, vcfs) {
		gr <- vcfs[[id]]
		if (!is.na((metadata %>% filter(Id == id))$CX_CALLER)) {
			# only looking at intrachromosomal calls
			gr <- gr[!is.na(gr$svLen),]
			#if (str_detect((metadata %>% filter(Id == id))$CX_REFERENCE_VCF_VARIANTS, "BP")) {
				# filtering small events calls to remove spurious indels caused by sequence homology around breakpoints
			#	gr <- gr[abs(gr$svLen) >= 50,]
			#}
		}
		return(gr)
	}

