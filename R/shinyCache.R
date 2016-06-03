source("sv_benchmark.R")
library(R.cache)

#' @param output of previous call to LoadPlotData. This will be reused as much as possible to avoid recalculation
LoadPlotData <- function(
		datadir,
		maxgap,
		ignore.strand,
		sizemargin,
		ignore.duplicates,
		ignore.interchromosomal,
		mineventsize,
		maxeventsize,
		vcftransform,
		truthgr,
		existingCache) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	# set up all cache keys for all the data
	truthhash <- NULL
	if (!is.null(truthgr)) {
	  truthhash <- getChecksum(truthgr)
	}
	keymetadata <- list(datadir)
	keyvcfs <- list(keymetadata)
	keycalls <- list(keyvcfs, maxgap, ignore.strand, sizemargin, truthhash) # TODO: vcftransform
	keydfs <- list(keycalls, ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize)
	slice <- list(
		datadir=datadir,
		keymetadata=keymetadata,
		keyvcfs=keyvcfs,
		keycalls=keycalls,
		keydfs=keydfs
		)
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
		slice$vcfs <- existingCache$vcfs
	}
	if (getChecksum(keycalls) == getChecksum(existingCache$keycalls)) {
		slice$calls <- existingCache$calls
	}
	if (getChecksum(keydfs) == getChecksum(existingCache$keydfs)) {
		slice$dfs <- existingCache$dfs
	}
	# Load metadata
	if (is.null(slice$metadata)) {
		slice$metadata <- loadCache(key=keymetadata, dirs=".Rcache/LoadPlotData/metadata")
		if (is.null(slice$metadata)) {
			slice$metadata <- LoadMetadata(datadir)
			saveCache(slice$metadata, key=keymetadata, dirs=".Rcache/LoadPlotData/metadata")
		}
	}
	# Load graphs, loading from R.cache whenever possible to avoid recalculation
	if (is.null(slice$dfs)) {
		slice$dfs <- loadCache(key=keydfs, dirs=".Rcache/LoadPlotData/dfs")
		write("Recalculating dfs", stderr())
		if (is.null(slice$dfs)) {
			# To recalculate the graph we need the call set
			if (is.null(slice$calls)) {
				# Load call set
				slice$calls <- loadCache(key=keycalls, dirs=".Rcache/LoadPlotData/calls")
				if (is.null(slice$calls)) {
					# To recalculate the call set we need the vcfs
					write("Recalculating calls", stderr())
					if (is.null(slice$vcfs)) {
						slice$vcfs <- loadCache(key=keyvcfs, dirs=".Rcache/LoadPlotData/vcfs")
						if (is.null(slice$vcfs)) {
							# metadata is always loaded so we're fine
							write("Loading vcfs", stderr())
							slice$vcfs <- LoadMinimalSVFromVCF(datadir, metadata=slice$metadata)
							saveCache(slice$vcfs, key=keyvcfs, dirs=".Rcache/LoadPlotData/vcfs")
						}
					}
					slice$calls <- LoadCallSets(slice$metadata, slice$vcfs, maxgap, ignore.strand, sizemargin, vcftransform, truthgr)
					saveCache(slice$calls, key=keycalls, dirs=".Rcache/LoadPlotData/calls")
				}
			}
			slice$dfs <- LoadGraphDataFrames(slice$metadata, slice$calls, ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize)
			saveCache(slice$dfs, key=keydfs, dirs=".Rcache/LoadPlotData/dfs")
		}
	}
	setCacheRootPath(cacheroot)
	return(slice)
}
LoadVCFs <- function(datadir, metadata) {
	vcfs <- LoadMinimalSVFromVCF(datadir, metadata=metadata)
	return(vcfs)
}
LoadCallSets <- function(metadata, vcfs, maxgap, ignore.strand, sizemargin, vcftransform, truthgr) {
	vcfs <- sapply(names(vcfs), vcftransform, metadata, vcfs, simplify=FALSE, USE.NAMES=TRUE)
	mcalls <- NULL
	for (includeFiltered in c(TRUE, FALSE)) {
		calls <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=includeFiltered, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, truthgr=truthgr)
		calls$truth <- calls$truth %>% mutate(duptp=FALSE)
		calls$calls <- calls$calls%>% filter(!tp) # take the truth vcf version of tp calls since it has the actual event size
		mergedcalls <- rbind(calls$calls, calls$truth)
		if (includeFiltered) {
			mergedcalls$CallSet <- "High & Low confidence"
		} else {
			mergedcalls$CallSet <- "High confidence only"
		}
		mcalls <- rbind(mcalls, mergedcalls)
	}
	return(mcalls)
}
LoadGraphDataFrames <- function(metadata, calls, ignore.duplicates, ignore.interchromosomal=TRUE, mineventsize=NULL, maxeventsize=NULL) {
	if (ignore.duplicates) {
		calls <- calls %>% filter(!duptp)
	}
	if (ignore.interchromosomal) {
		calls <- calls %>% filter(!is.na(svLen))
	}
	if (!is.null(mineventsize)) {
		calls <- calls %>% filter(svLen >= mineventsize)
	}
	if (!is.null(maxeventsize)) {
		calls <- calls %>% filter(svLen <= maxeventsize)
	}
	mostSensitiveAligner <- calls %>%
		select(Id, CallSet, tp) %>%
		group_by(Id, CallSet) %>%
		summarise(tp=sum(tp)) %>%
		ungroup() %>%
		arrange(desc(tp)) %>%
		left_join(metadata) %>%
		distinct(CallSet, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF) %>%
		select(Id, CallSet)
	callsByEventSize <- calls %>%
	  #filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
	  filter(!fp) %>%
	  filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	  group_by(Id, CallSet, svLen) %>%
	  summarise(sens=sum(tp)/sum(tp+fn)) %>%
	  ungroup() %>%
	  left_join(metadata) %>%
	  mutate(eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS))
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
	  left_join(metadata) %>%
	  mutate(eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS))
	return(list(mostSensitiveAligner=mostSensitiveAligner, callsByEventSize=callsByEventSize, roc=roc))
}
