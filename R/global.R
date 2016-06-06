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
