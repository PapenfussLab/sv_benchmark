# per-VCF version of shinyCache that minimises the R memory footprint by
# caching all information except the final on a per-VCF basis
source("sv_benchmark.R")
library(R.cache)
library(dplyr)

#'
LoadCachedMetadata <- function(datadir) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	keymetadata <- list(datadir)
	metadata <- loadCache(key=keymetadata, dirs=".Rcache/metadata")
	if (is.null(metadata)) {
		metadata <- LoadMetadata(datadir)
		saveCache(metadata, key=keymetadata, dirs=".Rcache/metadata")
	}
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
		existingCache,
		loadFromCacheOnly=TRUE,
		loadAll=FALSE) {
	if (!is.null(eventtypes)) {
		# order does not matter; strip names that were confusing the cache
		eventtypes <- as.character(sort(eventtypes))
	}
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)

	args$LoadCachedMetadata <- list(datadir=datadir)
	args$LoadLongReadTruthgr <- list(bedpedir=truthbedpedir, mintruthbedpescore=mintruthbedpescore)
	args$TransformVcf <- list(grtransform=grtransform, grtransformName=grtransformName)
	args$LoadCalls <- list(maxgap, ignore.strand, sizemargin, requiredHits)
	args$LoadGraphDataFrame <- list(ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes)

	metadata <- LoadCachedMetadata(datadir)
	slice <- list(
		dfs <- LoadCachedGraphDataFrame(args)
	)
	setCacheRootPath(cacheroot)
	return(slice)
}
LoadGraphDataFrame <- function(args) {
	# todo: cache
	# todo
	metadata <- LoadCachedMetadata(datadir)
	return(metadata %>%
		group_by(Id) %>%
		do(LoadGraphDataFrameForId(.$Id, args)))
}
LoadGraphDataFrameForId <- function(args, id, ignore.duplicates, ignore.interchromosomal=TRUE, mineventsize=NULL, maxeventsize=NULL, eventtypes=NULL, rocSlicePoints=100) {
	metadata <- LoadCachedMetadata(datadir)
	calls <- LoadCallsForId(id, args)

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
LoadCallsForId <- function(id, ) {
	# todo: remove all ScoreVariantsFromTruth caching
	# todo: load truthgr from metadata if truth arg is null
	callgr <- TransformVcf(args, id)
	if (!is.null(grtransform)) {
		callgrlist <- sapply(names(callgrlist), function(id, metadata, callgrlist) {
				return(grtransform(callgrlist[[id]], metadata %>% filter(Id == id)))
			}, metadata, callgrlist, simplify=FALSE, USE.NAMES=TRUE)
	}
	mcalls <- NULL
	for (includeFiltered in c(TRUE, FALSE)) {
		calls <- ScoreVariantsFromTruth(callgrlist, metadata, includeFiltered=includeFiltered, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgr=truthgr, keytruth=keytruth, keycalls=list(grtransformName))
		mergedcalls <- calls$calls
		if (!is.null(calls$truth)) {
			calls$truth <- calls$truth %>% mutate(duptp = FALSE)
			calls$calls <- calls$calls %>% filter(!tp) # take the truth vcf version of tp calls since it has the actual event size
			mergedcalls <- rbind(calls$calls, calls$truth)
		}
		mcalls <- rbind(mcalls, mergedcalls %>% mutate(CallSet = ifelse(includeFiltered, "High & Low confidence", "High confidence only")))
	}
	return(mcalls)
}
CachedTransformVcf <- function(id, grtransform, grtransformName) {
	return (TransformVcf(id, grtransformName))
}
TransformVcf <- function(id, grtransform) {
	list(filename, caller, transform, nominalPosition)
	vcf <- readVcf(filename, "hg19")
	if (!is.null(transform)) {
		vcf <- transform(vcf)
	}
	vcf <- withqual(vcf, caller)
	gr <- breakpointRanges(vcf, nominalPosition)
	gr$paramRangeID <- NULL
	gr$REF <- NULL
	gr$ALT <- NULL
	#gr$svtype <- NULL
	#gr$svLen <- NULL
	gr$insSeq <- NULL
	#gr$insLen <- NULL
	return(gr)
}

