source("sv_benchmark.R")
library(R.cache)

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
#' @param output of previous call to LoadPlotData. This will be reused as much as possible to avoid recalculation
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
		grtransform,
		truthgr,
		eventtypes,
		existingCache) {
  if (!is.null(eventtypes)) {
    # order does not matter; strip names that were confusing the cache
    eventtypes <- as.character(sort(eventtypes))
  }
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	# set up all cache keys for all the data
	truthhash <- NULL
	if (!is.null(truthgr)) {
		truthhash <- getChecksum(truthgr)
	}
	keymetadata <- list(datadir)
	keyvcfs <- list(keymetadata)
	keycalls <- list(keyvcfs, maxgap, ignore.strand, sizemargin, requiredHits, truthhash, grtransform)
	keydfs <- list(keycalls, ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes)
	slice <- list(
		datadir=datadir,
		keymetadata=keymetadata,
		keyvcfs=keyvcfs,
		keycalls=keycalls,
		keydfs=keydfs
		)
	# # Debug cache missing
	# write(sprintf("datadir=%s", datadir), stderr())
	# write(sprintf("keymetadata=%s", getChecksum(keymetadata)), stderr())
	# write(sprintf("keyvcfs=%s", getChecksum(keyvcfs)), stderr())
	# write(sprintf("keycalls=%s", getChecksum(keycalls)), stderr())
	# write(sprintf(".keyvcfs=%s", getChecksum(keyvcfs)), stderr())
	# write(sprintf(".maxgap=%s", getChecksum(maxgap)), stderr())
	# write(sprintf(".ignore.strand=%s", getChecksum(ignore.strand)), stderr())
	# write(sprintf(".sizemargin=%s", getChecksum(sizemargin)), stderr())
	# write(sprintf(".requiredHits=%s:%s", getChecksum(requiredHits), requiredHits), stderr())
	# write(sprintf(".truthhash=%s", getChecksum(truthhash)), stderr())
	# write(sprintf(".grtransform=%s", getChecksum(grtransform)), stderr())
	# write(sprintf("keydfs=%s", getChecksum(keydfs)), stderr())
	# stop()
	if (is.null(existingCache)) {
		existingCache <- slice
	}
	# discard cache if using different data
	if (existingCache$datadir != slice$datadir) {
		existingCache <- slice
	}
	# Copy in-memory cached objects that are still valid
	if (getChecksum(keymetadata) == getChecksum(existingCache$keymetadata)) {
		slice$metadata <- existingCache$metadata
	}
	if (getChecksum(keyvcfs) == getChecksum(existingCache$keyvcfs)) {
		slice$callgrlist <- existingCache$callgrlist
	}
	if (getChecksum(keycalls) == getChecksum(existingCache$keycalls)) {
		slice$calls <- existingCache$calls
	}
	if (getChecksum(keydfs) == getChecksum(existingCache$keydfs)) {
		slice$dfs <- existingCache$dfs
	}
	# Load metadata
	if (is.null(slice$metadata)) {
		slice$metadata <- LoadCachedMetadata(datadir)
	}
	# Load graphs, loading from R.cache whenever possible to avoid recalculation
	if (is.null(slice$dfs)) {
	  write(sprintf("Loading dfs (%s)", getChecksum(keydfs)), stderr())
		slice$dfs <- loadCache(key=keydfs, dirs=".Rcache/LoadPlotData/dfs")
		if (is.null(slice$dfs)) {
			write(sprintf("Recalculating dfs (%s)", getChecksum(keydfs)), stderr())
			# To recalculate the graph we need the call set
			if (is.null(slice$calls)) {
				# Load call set
				slice$calls <- loadCache(key=keycalls, dirs=".Rcache/LoadPlotData/calldf")
				if (is.null(slice$calls)) {
					write(sprintf("Recalculating calls (%s)", getChecksum(keycalls)), stderr())
					# To recalculate the call set we need the SV grs from the vcfs
					if (is.null(slice$callgrlist)) {
						slice$callgrlist <- loadCache(key=keyvcfs, dirs=".Rcache/LoadPlotData/callgr")
						if (is.null(slice$callgrlist)) {
							# metadata is always loaded so we're fine
						  write(sprintf("Loading vcfs (%s)", getChecksum(keyvcfs)), stderr())
							slice$callgrlist <- LoadVCFs(datadir, metadata=slice$metadata)
							saveCache(slice$callgrlist, key=keyvcfs, dirs=".Rcache/LoadPlotData/callgr")
						}
					}
					slice$calls <- LoadCallSets(slice$metadata, slice$callgrlist, maxgap, ignore.strand, sizemargin, requiredHits, grtransform, truthgr)
					saveCache(slice$calls, key=keycalls, dirs=".Rcache/LoadPlotData/calldf")
				}
			}
			slice$dfs <- LoadGraphDataFrames(slice$metadata, slice$calls, ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes)
			saveCache(slice$dfs, key=keydfs, dirs=".Rcache/LoadPlotData/dfs")
		}
	}
	setCacheRootPath(cacheroot)
	return(slice)
}
LoadVCFs <- function(datadir, metadata) {
	callgrlist <- LoadMinimalSVFromVCF(datadir, metadata=metadata)
	return(callgrlist)
}
LoadCallSets <- function(metadata, callgrlist, maxgap, ignore.strand, sizemargin, requiredHits, grtransform, truthgr) {
	if (!is.null(grtransform)) {
		callgrlist <- sapply(names(callgrlist), function(id, metadata, callgrlist) {
				return(grtransform(callgrlist[[id]], metadata %>% filter(Id == id)))
			}, metadata, callgrlist, simplify=FALSE, USE.NAMES=TRUE)
	}
	mcalls <- NULL
	for (includeFiltered in c(TRUE, FALSE)) {
		calls <- ScoreVariantsFromTruth(callgrlist, metadata, includeFiltered=includeFiltered, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgr=truthgr)
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
LoadGraphDataFrames <- function(metadata, calls, ignore.duplicates, ignore.interchromosomal=TRUE, mineventsize=NULL, maxeventsize=NULL, eventtypes=NULL) {
	if (ignore.duplicates) {
		calls <- calls %>% filter(!duptp)
	}
	if (ignore.interchromosomal) {
		calls <- calls %>% filter(!is.na(svLen))
	}
	if (!is.null(mineventsize)) {
		calls <- calls %>% filter(abs(svLen + insLen) >= mineventsize)
	}
	if (!is.null(maxeventsize)) {
		calls <- calls %>% filter(abs(svLen + insLen) <= maxeventsize)
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
		left_join(md) %>%
		distinct(CallSet, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF, .keep_all = TRUE) %>%
		select(Id, CallSet)
	callsByEventSize <- calls %>%
		#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
		filter(!fp) %>%
		#filter(Id %in% md$Id[md$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
		group_by(Id, CallSet, svLen) %>%
		summarise(sens=sum(tp)/sum(tp+fn)) %>%
		ungroup() %>%
		left_join(md)
	roc <- calls %>%
		select(Id, CallSet, QUAL, tp, fp, fn) %>%
		rbind(calls %>% select(Id, CallSet) %>% distinct(Id, CallSet) %>% mutate(QUAL=max(calls$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
		#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
		group_by(Id, CallSet) %>%
		arrange(desc(QUAL)) %>%
		mutate(events=sum(tp) + sum(fn)) %>%
		mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
		group_by(Id, CallSet, QUAL, events) %>%
		summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
		ungroup() %>%
		mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
		left_join(md)
	# lossless reduction of roc plot points by elimination of points on straight lines
	roc <- roc %>%
		group_by(Id, CallSet) %>%
		arrange(tp + fp) %>%
		filter(
			# keep start/end
			is.na(lag(tp)) | is.na(lead(tp)) |
			# keep group transitions (TODO: is there a way to make lead/lag across group_by return NA?)
			Id != lag(Id) | CallSet != lag(CallSet) |
			Id != lead(Id) | CallSet != lead(CallSet) |
			# slopes not equal dx1/dy1 != dx2/dy2 -> dx1*dy2 != dx2*dy1
			(tp - lag(tp))*(lead(fp) - lag(fp)) != (lead(tp) - lag(tp))*(fp - lag(fp))) %>%
		ungroup()
	# lossy removal of points with least change
	for (n in c(4, 16, 32, 64)) {
		roc <- roc %>%
			group_by(Id, CallSet) %>%
			arrange(tp + fp) %>%
			filter(
				is.na(lag(tp)) | is.na(lead(tp)) |
				Id != lag(Id) | CallSet != lag(CallSet) |
				Id != lead(Id) | CallSet != lead(CallSet) |
				# remove points with least amount of change
				lead(tp) - lag(tp) + lead(fp) - lag(fp) > n |
				# keep every 5th to prevent removal of large segments
				row_number() %% 5 == 0
			) %>%
			ungroup()
	}
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
		ungroup() %>%
		mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
		left_join(md)
	# lossless reduction of roc plot points by elimination of points on straight lines
	rocbyrepeat <- rocbyrepeat %>%
		group_by(Id, CallSet, repeatClass) %>%
		arrange(tp + fp) %>%
		filter(
			# keep start/end
			is.na(lag(tp)) | is.na(lead(tp)) |
			# keep group transitions (TODO: is there a way to make lead/lag across group_by return NA?)
			Id != lag(Id) | CallSet != lag(CallSet) | repeatClass != lag(repeatClass) |
			Id != lead(Id) | CallSet != lead(CallSet) | repeatClass != lead(repeatClass) |
			# slopes not equal dx1/dy1 != dx2/dy2 -> dx1*dy2 != dx2*dy1
			(tp - lag(tp))*(lead(fp) - lag(fp)) != (lead(tp) - lag(tp))*(fp - lag(fp))) %>%
		ungroup()
	# lossy removal of points with least change
	for (n in c(4, 16, 32, 64)) {
		rocbyrepeat <- rocbyrepeat %>%
			group_by(Id, CallSet, repeatClass) %>%
			arrange(tp + fp) %>%
			filter(
				is.na(lag(tp)) | is.na(lead(tp)) |
				Id != lag(Id) | CallSet != lag(CallSet) | repeatClass != lag(repeatClass) |
				Id != lead(Id) | CallSet != lead(CallSet) | repeatClass != lead(repeatClass) |
				# remove points with least amount of change
				lead(tp) - lag(tp) + lead(fp) - lag(fp) > n |
				# keep every 5th to prevent removal of large segments
				row_number() %% 5 == 0
			) %>%
			ungroup()
	}

	bpErrorDistribution <- calls %>%
		filter(tp) %>%
		select(Id, CallSet, bperror) %>%
		group_by(Id, CallSet, bperror) %>%
		summarize(n=n())
	bpErrorDistribution <- bpErrorDistribution %>%
		left_join(bpErrorDistribution %>% group_by(Id, CallSet) %>% summarize(count=sum(n))) %>%
		mutate(rate=n/count) %>%
		select(-count) %>%
		left_join(md)
	return(list(mostSensitiveAligner=mostSensitiveAligner,
							callsByEventSize=callsByEventSize,
							roc=roc,
							rocbyrepeat=rocbyrepeat,
							bpErrorDistribution=bpErrorDistribution))
}
LoadLongReadTruthgr <- function(dir) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(dir)
	key <- list(dir=dir)
	gr <- loadCache(key = key, dirs = ".Rcache/LoadLongReadTruthgr")
	if (is.null(gr)) {
		gr <- .LoadLongReadTruthgr(dir)
		saveCache(gr, key = key, dirs = ".Rcache/LoadLongReadTruthgr")
	}
	setCacheRootPath(cacheroot)
	return(gr)
}
.LoadLongReadTruthgr <- function(dir) {
	gr <- NULL
	for (file in list.files(path = dir, pattern = ".bedpe.gz", full.names = TRUE)) {
		gr2 <- import.sv.bedpe(file)
		if (is.null(gr)) {
			gr <- gr2
		} else {
			gr <- c(gr, gr2)
		}
	}
	seqlevelsStyle(gr) <- "UCSC"
	return (gr)
}
