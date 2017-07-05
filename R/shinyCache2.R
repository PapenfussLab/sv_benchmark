# per-VCF version of shinyCache that minimises the R memory footprint by
# caching all information except the final on a per-VCF basis
source("sv_benchmark.R")
library(R.cache)
library(dplyr)
library(assertthat)

memcache.metadata <- list(key=NULL, value=NULL)
LoadCachedMetadata <- function(datadir) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	keymetadata <- list(datadir)
	# try in-memory cache
	if (getChecksum(memcache.metadata$key) == getChecksum(keymetadata)) {
		metadata <- memcache.metadata$value
		return (metadata)
	}
	metadata <- loadCache(key=keymetadata, dirs=".Rcache/metadata")
	if (is.null(metadata)) {
		write(sprintf("LoadMetadata %s (%s)", datadir, getChecksum(datadir)), stderr())
		metadata <- LoadMetadata(datadir)
		if (is.null(metadata$CX_REFERENCE_VCF)) {
			# Hack to force a default truth even if none exists
			metadata$CX_REFERENCE_VCF <- "00000000000000000000000000000002.reference.vcf"
		}
		saveCache(metadata, key=keymetadata, dirs=".Rcache/metadata")
	}
	memcache.metadata$key <<- keymetadata
	memcache.metadata$value <<- metadata
	setCacheRootPath(cacheroot)
	return(metadata)
}
LoadPlotData <- function(
		datadir,
		maxgap,
		ignore.strand,
		sizemargin,
		ignore.duplicates,
		ignore.interchromosomal,
		requiredHits,
		mineventsize,
		maxeventsize,
		grtransformName,
		grtransform,
		truthbedpedir,
		mintruthbedpescore,
		eventtypes,
		# ignored: only used for compatability with shinyCache.R implementation
		existingCache,
		loadFromCacheOnly=TRUE,
		loadAll=FALSE) {
	if (!is.null(eventtypes)) {
		# order does not matter; strip names that were confusing the cache
		eventtypes <- as.character(sort(eventtypes))
	}
	# ensure numeric values are of the correct type
	if (!is.null(maxgap)) {
		maxgap <- as.numeric(maxgap)
	}
	if (!is.null(sizemargin)) {
		sizemargin <- as.numeric(sizemargin)
	}
	if (!is.null(requiredHits)) {
		requiredHits <- as.numeric(requiredHits)
	}
	if (!is.null(mineventsize)) {
		mineventsize <- as.numeric(mineventsize)
	}
	if (!is.null(maxeventsize)) {
		maxeventsize <- as.numeric(maxeventsize)
	}
	if (!is.null(mintruthbedpescore)) {
		mintruthbedpescore <- as.numeric(mintruthbedpescore)
	}
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	metadata <- LoadCachedMetadata(datadir)
	truthgr <- NULL
	truthgrName <- ""
	if (!is.null(truthbedpedir)) {
		truthgrName <- truthbedpedir
	}
	rocSlicePoints <- 100
	nominalPosition <- FALSE
	slice <-list()
	slice$metadata <- metadata
	slice$dfs <- .CachedLoadGraphDataFrame(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata,
		maxgap, sizemargin, ignore.strand, requiredHits,
			# use lazy evaluation for truthgr and truthgrName
			if (is.null(truthbedpedir)) { NULL } else {.CachedLoadTruthBedpe(truthbedpedir, mintruthbedpescore)},
			if (is.null(truthbedpedir)) { NULL } else {truthbedpedir},
		grtransform, grtransformName, nominalPosition)

	setCacheRootPath(cacheroot)
	return(slice)
}
.CachedLoadTruthBedpe <- function(truthbedpedir, mintruthbedpescore) {
	cachekey <- list(truthbedpedir, mintruthbedpescore)
	cachedir <- ".Rcache/TruthBedpe"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadTruthBedpe %s (%s)", truthbedpedir, getChecksum(cachekey)), stderr())
		result <- .LoadTruthBedpe(truthbedpedir, mintruthbedpescore)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
import.sv.bedpe.dir <- function(dir) {
	gr <- NULL
	for (file in list.files(path = dir, pattern = ".bedpe.gz", full.names = TRUE)) {
		gr2 <- import.sv.bedpe(file)
		if (is.null(gr)) {
			gr <- gr2
		} else {
			gr <- c(gr, gr2)
		}
	}
	return (gr)
}
.LoadTruthBedpe <- function(truthbedpedir, mintruthbedpescore) {
	truthgr <- import.sv.bedpe.dir(truthbedpedir)
	truthgr <- truthgr[truthgr$score >= mintruthbedpescore]
	seqlevelsStyle(truthgr) <- "UCSC"
	return(truthgr)
}
.CachedLoadGraphDataFrame <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {
	cachekey <- cachekey <- list(
		ignore.duplicates=ignore.duplicates, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, eventtypes=eventtypes, rocSlicePoints=rocSlicePoints,
		datadir=datadir,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgrName=truthgrName,
		grtransformName=grtransformName, nominalPosition=nominalPosition)
	cachedir <- ".Rcache/GraphDataFrame"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadGraphDataFrame %s (%s)", datadir, getChecksum(cachekey)), stderr())
		result <- .LoadGraphDataFrame(
			ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
			datadir, metadata,
			maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
			grtransform, grtransformName, nominalPosition)
		saveCache(result, key=cachekey, dirs=cachedir)
	}
	return(result)
}
.LoadGraphDataFrame <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {
	ids <- (metadata %>% filter(!is.na(CX_CALLER)))$Id
	ids <- sample(ids) # randomise our processing order so parallel caches won't all try to load the same record at the same time
	dfslist <- lapply(ids, function(id)
		.CachedLoadGraphDataFrameForId(
			ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
			datadir, metadata, id,
			maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
			grtransform, grtransformName, nominalPosition)
		)
	return(list(
		mostSensitiveAligner=bind_rows(lapply(dfslist, function(item) item$mostSensitiveAligner)),
		callsByEventSize=bind_rows(lapply(dfslist, function(item) item$callsByEventSize)),
		roc=bind_rows(lapply(dfslist, function(item) item$roc)),
		rocbyrepeat=bind_rows(lapply(dfslist, function(item) item$rocbyrepeat)),
		bpErrorDistribution=bind_rows(lapply(dfslist, function(item) item$bpErrorDistribution))
	))
}
.CachedLoadGraphDataFrameForId <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {
	cachekey <- list(
		ignore.duplicates=ignore.duplicates, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, eventtypes=eventtypes, rocSlicePoints=rocSlicePoints,
		datadir=datadir, id=id,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgrName=truthgrName,
		grtransformName=grtransformName, nominalPosition=nominalPosition)
	cachedir <- ".Rcache/GraphDataFrameById"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadGraphDataFrameForId %s (%s)", id, getChecksum(cachekey)), stderr())
		nullwrap <- function(x) { if (is.null(x)) { return("NULL") } else { return(x) } }
		write(paste0(".LoadGraphDataFrameForId(",
								 "ignore.duplicates=", ignore.duplicates, ",",
								 "ignore.interchromosomal=", ignore.interchromosomal, ",",
								 "mineventsize=", nullwrap(mineventsize), ",",
								 "maxeventsize=", nullwrap(maxeventsize), ",",
								 "eventtypes=", nullwrap(eventtypes), ",",
								 "rocSlicePoints=", rocSlicePoints, ",",
								 "datadir=\"", datadir, "\",",
								 "metadata=LoadCachedMetadata(\"", datadir, "\"),",
								 "id=\"", id, "\",",
								 "maxgap=", maxgap, ",",
								 "sizemargin=", sizemargin, ",",
								 "ignore.strand=", ignore.strand, ",",
								 "requiredHits=", requiredHits, ",",
								 "truthgr=truthgr,",
								 "truthgrName=\"", truthgrName, "\",",
								 "grtransform=simoptions$grtransform[[\"",grtransformName,"\"]],",
								 "grtransformName=\"", grtransformName, "\",",
								 "nominalPosition=", nominalPosition, ")"
								 ), stderr())
		result <- .LoadGraphDataFrameForId(
			ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
			datadir, metadata, id,
			maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
			grtransform, grtransformName, nominalPosition)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.LoadGraphDataFrameForId <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {

	calls <- .CachedLoadCallsForId(datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition)
	if (is.null(calls)) {
		return(NULL)
	}
	if (ignore.duplicates) {
		calls <- calls %>% filter(!duptp)
	}
	if (ignore.interchromosomal) {
		calls <- calls %>% filter(!is.na(svLen))
	}
	if (!is.null(mineventsize)) {
		calls <- calls %>% filter(is.na(svLen) | abs(svLen + insLen) >= mineventsize)
	}
	if (!is.null(maxeventsize)) {
		calls <- calls %>% filter(is.na(svLen) | abs(svLen + insLen) <= maxeventsize)
	}
	if (is.null(metadata$CX_MULTIMAPPING_LOCATIONS)) {
		metadata$CX_MULTIMAPPING_LOCATIONS <- NA_integer_
	}
	if (!is.null(eventtypes)) {
		calls <- calls %>% filter(simpleEvent %in% eventtypes)
	}

	md <- metadata %>%
		select(Id, CX_ALIGNER, CX_ALIGNER_MODE, CX_MULTIMAPPING_LOCATIONS, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH)
	if (!is.null(metadata$CX_REFERENCE_VCF_VARIANTS)) {
		md$CX_REFERENCE_VCF_VARIANTS <- metadata$CX_REFERENCE_VCF_VARIANTS
		md$eventtype <- PrettyVariants(md$CX_REFERENCE_VCF_VARIANTS)
	}
	if (!is.null(metadata$CX_REFERENCE_VCF)) {
		md$CX_REFERENCE_VCF <- metadata$CX_REFERENCE_VCF
	} else {
		md$CX_REFERENCE_VCF <- "longread"
	}
	mostSensitiveAligner <- calls %>%
		select(Id, CallSet, tp) %>%
		group_by(Id, CallSet) %>%
		summarise(tp=sum(tp)) %>%
		ungroup() %>%
		arrange(desc(tp)) %>%
		left_join(md, by="Id") %>%
		distinct(CallSet, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF, .keep_all = TRUE) %>%
		select(Id, CallSet)
	callsByEventSize <- NULL
	if (!is.null(md$eventtype)) { # only generate callsByEventSize for simulation data sets
		callsByEventSize <- calls %>%
			#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
			filter(!fp) %>%
			group_by(Id, CallSet, svLen) %>%
			summarise(sens=sum(tp)/sum(tp+fn)) %>%
			ungroup() %>%
			left_join(md) %>%
			# Don't include "Breakpoint","SINE/ALU Breakpoint" in the simulation
			filter(eventtype %in% c("Tandem Duplication","Deletion","Insertion","Inversion"))
	}
	roc <- calls %>%
		select(Id, CallSet, QUAL, tp, fp, fn) %>%
		rbind(calls %>% select(Id, CallSet) %>% distinct(Id, CallSet) %>% mutate(QUAL=max(calls$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
		#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
		group_by(Id, CallSet) %>%
		arrange(desc(QUAL)) %>%
		mutate(events=sum(tp) + sum(fn)) %>%
		mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
		# each QUAL score is a point on the ROC plott
		group_by(Id, CallSet, QUAL, events) %>%
		summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
		# QUAL scores with the same number of tp calls can be merged on the ROC plot
		group_by(Id, CallSet, tp, events) %>%
		summarise(fp=max(fp), fn=max(fn), QUAL=min(QUAL)) %>%
		# subsample along tp and tp+fp axis
		group_by(Id, CallSet) %>%
		slice(unique(c(
			1,
			findInterval(seq(0, max(tp), max(tp)/rocSlicePoints), tp),
			findInterval(seq(0, max(tp + fp), max(tp + fp)/rocSlicePoints), tp + fp),
			n()
		))) %>%
		ungroup() %>%
		mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
		left_join(md, by="Id")

	rocbyrepeat <- calls %>%
		select(Id, CallSet, repeatClass, QUAL, tp, fp, fn) %>%
		rbind(calls %>% select(Id, CallSet, repeatClass) %>% distinct(Id, CallSet, repeatClass) %>% mutate(QUAL=max(calls$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
		#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
		group_by(Id, CallSet, repeatClass) %>%
		arrange(desc(QUAL)) %>%
		mutate(events=sum(tp) + sum(fn)) %>%
		mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
		group_by(Id, CallSet, repeatClass, QUAL, events) %>%
		summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
		group_by(Id, CallSet, repeatClass, tp, events) %>%
		summarise(fp=max(fp), fn=max(fn), QUAL=min(QUAL)) %>%
		group_by(Id, CallSet, repeatClass) %>%
		# subsample along tp and tp+fp axis
		slice(unique(c(
			1,
			findInterval(seq(0, max(tp), max(tp)/rocSlicePoints), tp),
			findInterval(seq(0, max(tp + fp), max(tp + fp)/rocSlicePoints), tp + fp),
			n()
		))) %>%
		ungroup() %>%
		mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
		left_join(md, by="Id")
	# lossless reduction of roc plot points by elimination of points on straight lines
	# rocbyrepeat <- rocbyrepeat %>%
	# 	group_by(Id, CallSet, repeatClass) %>%
	# 	arrange(tp + fp) %>%
	# 	filter(
	# 		# keep start/end
	# 		is.na(lag(tp)) | is.na(lead(tp)) |
	# 		# keep group transitions (TODO: is there a way to make lead/lag across group_by return NA?)
	# 		Id != lag(Id) | CallSet != lag(CallSet) | repeatClass != lag(repeatClass) |
	# 		Id != lead(Id) | CallSet != lead(CallSet) | repeatClass != lead(repeatClass) |
	# 		# slopes not equal dx1/dy1 != dx2/dy2 -> dx1*dy2 != dx2*dy1
	# 		(tp - lag(tp))*(lead(fp) - lag(fp)) != (lead(tp) - lag(tp))*(fp - lag(fp))) %>%
	# 	ungroup()
	# lossy removal of points with least change
	# for (n in c(4, 16, 32, 64)) {
	# 	rocbyrepeat <- rocbyrepeat %>%
	# 		group_by(Id, CallSet, repeatClass) %>%
	# 		arrange(tp + fp) %>%
	# 		filter(
	# 			is.na(lag(tp)) | is.na(lead(tp)) |
	# 			Id != lag(Id) | CallSet != lag(CallSet) | repeatClass != lag(repeatClass) |
	# 			Id != lead(Id) | CallSet != lead(CallSet) | repeatClass != lead(repeatClass) |
	# 			# remove points with least amount of change
	# 			lead(tp) - lag(tp) + lead(fp) - lag(fp) > n |
	# 			# keep every 5th to prevent removal of large segments
	# 			row_number() %% 5 == 0
	# 		) %>%
	# 		ungroup()
	# }

	bpErrorDistribution <- calls %>%
		filter(tp) %>%
		select(Id, CallSet, bperror) %>%
		group_by(Id, CallSet, bperror) %>%
		summarize(n=n())
	bpErrorDistribution <- bpErrorDistribution %>%
		left_join(bpErrorDistribution %>% group_by(Id, CallSet) %>% summarize(count=sum(n))) %>%
		mutate(rate=n/count) %>%
		select(-count) %>%
		left_join(md, by="Id")
	return(list(mostSensitiveAligner=mostSensitiveAligner,
							callsByEventSize=callsByEventSize,
							roc=roc,
							rocbyrepeat=rocbyrepeat,
							bpErrorDistribution=bpErrorDistribution))
}
.CachedLoadCallsForId <- function(datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {
	cachekey <- list(
		datadir=datadir, id=id,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgrName=truthgrName,
		grtransformName=grtransformName, nominalPosition=nominalPosition)
	cachedir <- ".Rcache/Calls"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadCallsForId %s (%s)", id, getChecksum(cachekey)), stderr())
		result <- .LoadCallsForId(
			datadir=datadir, metadata=metadata, id=id,
			maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgr=truthgr,
			grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.LoadCallsForId <- function(datadir, metadata, id,
		# ScoreVariantsFromTruth
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr,
		# .CachedTransformVcf
		grtransform, grtransformName, nominalPosition) {
	if (is.null(truthgr)) {
		truthid <- GetId((metadata %>% filter(Id==id))$CX_REFERENCE_VCF)
		if (is.na(truthid) | is.null(truthid)) {
			stop("Missing truth for ", id)
		}
		truthgr <- .CachedTransformVcf(datadir=datadir, metadata=metadata, id=truthid, grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
	}
	if (is.null(truthgr)) {
		stop("Missing truth for ", id)
	}
	callgr <- .CachedTransformVcf(datadir=datadir, metadata=metadata, id=id, grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
	if (is.null(callgr)) {
		return(NULL)
	}
	mcalls <- NULL
	for (includeFiltered in c(TRUE, FALSE)) {
		calls <- .ScoreVariantsFromTruthVCF(callgr=callgr, truthgr=truthgr, includeFiltered=includeFiltered, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, id=id, requiredHits=requiredHits)
		mergedcalls <- calls$calls
		if (!is.null(calls$truth)) {
			calls$truth <- calls$truth %>% mutate(duptp=FALSE)
			calls$calls <- calls$calls %>% filter(!tp) # take the truth vcf version of tp calls since it has the actual event size
			mergedcalls <- rbind(calls$calls, calls$truth)
		}
		mcalls <- rbind(mcalls, mergedcalls %>% mutate(CallSet = ifelse(includeFiltered, "High & Low confidence", "High confidence only")))
	}
	return(mcalls)
}
.CachedTransformVcf <- function(datadir, metadata, id, grtransform, grtransformName, nominalPosition) {
	cachekey <- list(datadir, id, grtransformName, nominalPosition)
	cachedir <- ".Rcache/vcfgr"
	assert_that(!is.null(id))
	assert_that(!is.na(id))
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		md <- metadata %>% filter(Id == id)
		if (nrow(md) != 1) {
			stop(paste(id, "not found in metadata"))
		}
		caller <- md$CX_CALLER
		write(sprintf(".TransformVcf %s (%s)", id, getChecksum(cachekey)), stderr())
		result <- .TransformVcf(datadir=datadir, id=id, caller=caller, grtransform=grtransform, nominalPosition=nominalPosition)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.TransformVcf <- function(datadir, id, caller, grtransform=NULL, vcftransform=NULL, nominalPosition) {
	#' Loads structural variant GRanges from the VCFs in the given directory
	filename <- list.files(datadir, pattern=paste0("^", id, ".*.vcf$"), full.names=TRUE)
	if (length(filename) == 0) {
		warning(paste("No vcf for id:", id, "skipping"))
		return(NULL)
	}
	if (length(filename) > 1) {
		stop(paste0("Multiple VCFs found for ", id))
	}
	vcf <- readVcf(filename, "")
	if (!is.null(vcftransform)) {
		vcf <- vcftransform(vcf)
	}
	vcf <- withqual(vcf, caller)
	gr <- breakpointRanges(vcf, nominalPosition)
	if (!is.null(grtransform)) {
		gr <- grtransform(gr)
	}
	gr$paramRangeID <- NULL
	gr$REF <- NULL
	gr$ALT <- NULL
	#gr$svtype <- NULL
	#gr$svLen <- NULL
	gr$insSeq <- NULL
	#gr$insLen <- NULL
	return(gr)
}

