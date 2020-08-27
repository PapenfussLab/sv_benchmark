# per-VCF version of shinyCache that minimises the R memory footprint by
# caching all information except the final on a per-VCF basis
source("sv_benchmark.R")
library(stringr)
library(R.cache)
library(dplyr)
library(data.table)
library(assertthat)

memcache.metadata <- list(key=NULL, value=NULL)
LoadCachedMetadata <- function(datadir) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	keymetadata <- list(datadir, files=list.files(datadir, pattern="*.metadata"))
	# try in-memory cache
	if (getChecksum(memcache.metadata$key) == getChecksum(keymetadata)) {
		metadata <- memcache.metadata$value
		return (metadata)
	}
	metadata <- loadCache(key=keymetadata, dirs=".Rcache/metadata")
	if (is.null(metadata)) {
		write(sprintf("LoadMetadata %s (%s)", datadir, getChecksum(datadir)), stderr())
		metadata <- LoadMetadata(datadir)
		if (str_detect(datadir, "na12878$")) {
			if (is.null(metadata$CX_REFERENCE_VCF)) {
				# Hack to force a default truth even if none exists
				metadata$CX_REFERENCE_VCF <- "00000000000000000000000000000002.reference.vcf" # TODO: why is the broken?
			}
			if (is.null(metadata$CX_SNP_TRUTH)) {
				# use the bwa/bcftools calls
				metadata$CX_SNP_TRUTH <- "2f4d15f5f6e428fb91e12dac971571b9.vcf"
			}
		} else if (str_detect(datadir, "chm13$")) {
			metadata$CX_REFERENCE_VCF <- "00000000000000000000000000000013.vcf"
			metadata$CX_SNP_TRUTH <- "c639766990fcfbca2eb45f7806362fe6.vcf"

		} else if (str_detect(datadir, "chm1$")) {
			metadata$CX_REFERENCE_VCF <- "00000000000000000000000000000001.vcf"
			metadata$CX_SNP_TRUTH <- "94bb6deef9f1bf1f9027a47e8488ae4f.vcf"

		} else if (str_detect(datadir, "chmboth$")) {
			metadata$CX_REFERENCE_VCF <- "00000000000000000000000000000014.vcf"
			metadata$CX_SNP_TRUTH <- "90f54397059a02d31adc4925ad57c439.vcf"

		} else if (str_detect(datadir, "HG002$")) {
			metadata$CX_REFERENCE_VCF <- "00000000000000000000000000000002.vcf"
			metadata$CX_SNP_TRUTH <- "2871cbd704c55d40fef7068e0a21dff9.vcf"
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
	eventtypes <- .cacheableEventTypes(eventtypes)
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
	truthgrName <- ""
	if (!is.null(truthbedpedir)) {
		truthgrName <- truthbedpedir
	}
	rocSlicePoints <- 100
	slice <-list()
	slice$metadata <- metadata
	slice$dfs <- .CachedLoadGraphDataFrame(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata,
		maxgap, sizemargin, ignore.strand, requiredHits,
			# use lazy evaluation for truthgr and truthgrName
			if (is.null(truthbedpedir)) { NULL } else {.CachedLoadTruthBedpe(truthbedpedir, mintruthbedpescore)},
			if (is.null(truthbedpedir)) { NULL } else {truthbedpedir},
		grtransform, grtransformName)

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
		grtransform, grtransformName) {
	cachekey <- list(
		ignore.duplicates=ignore.duplicates, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, eventtypes=eventtypes, rocSlicePoints=rocSlicePoints,
		datadir=datadir,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgrName=truthgrName,
		grtransformName=grtransformName)
	cachedir <- ".Rcache/GraphDataFrame"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadGraphDataFrame %s (%s)", datadir, getChecksum(cachekey)), stderr())
		result <- .LoadGraphDataFrame(
			ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
			datadir, metadata,
			maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
			grtransform, grtransformName)
		saveCache(result, key=cachekey, dirs=cachedir)
	}
	return(result)
}
.LoadGraphDataFrame <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName) {
	ids <- (metadata %>% filter(!is.na(CX_CALLER)))$Id
	ids <- sample(ids) # randomise our processing order so parallel caches won't all try to load the same record at the same time
	dfslist <- lapply(ids, function(id)
		.CachedLoadGraphDataFrameForId(
			ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
			datadir, metadata, id,
			maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
			grtransform, grtransformName)
		)
	# merge into single list
	result <- list()
	for (plotName in names(dfslist[[1]])) {
		result[[plotName]] <- bind_rows(lapply(dfslist, function(item) item[[plotName]]))
	}
	result$mostSensitiveAligner <- result$mostSensitiveAligner %>%
		arrange(desc(tp)) %>%
		distinct(CallSet, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF, .keep_all = TRUE) %>%
		dplyr::select(Id, CallSet)
	# Pair-wise correlation of callers?
	# Qual scores for each caller

	return(result)
}
.CachedLoadGraphDataFrameForId <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
		datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName) {
	cachekey <- list(
		ignore.duplicates=ignore.duplicates, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, eventtypes=eventtypes, rocSlicePoints=rocSlicePoints,
		datadir=datadir, id=id,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgrName=truthgrName,
		grtransformName=grtransformName)
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
								 "grtransformName=\"", grtransformName, "\")"
								 ), stderr())
		result <- .LoadGraphDataFrameForId(
			ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes, rocSlicePoints,
			datadir, metadata, id,
			maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
			grtransform, grtransformName)
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
		grtransform, grtransformName) {
	# browser()
	calls <- .LoadMinCallsForId(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes,
		datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, FALSE)
	nominalCalls <- .LoadMinCallsForId(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes,
		datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, TRUE)
	if (is.null(calls)) {
		return(NULL)
	}
	calls <- calls %>% mutate(nominalPosition=FALSE)
	nominalCalls <- nominalCalls %>% mutate(nominalPosition=FALSE)
	md <- metadata %>%
		dplyr::select(Id, CX_ALIGNER, CX_ALIGNER_MODE, CX_MULTIMAPPING_LOCATIONS, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH)
	if (!is.null(metadata$CX_REFERENCE_VCF_VARIANTS)) {
		md$CX_REFERENCE_VCF_VARIANTS <- metadata$CX_REFERENCE_VCF_VARIANTS
		md$eventtype <- PrettyVariants(md$CX_REFERENCE_VCF_VARIANTS)
	}
	if (!is.null(metadata$CX_REFERENCE_VCF)) {
		md$CX_REFERENCE_VCF <- metadata$CX_REFERENCE_VCF
	} else {
		md$CX_REFERENCE_VCF <- "longread"
	}
	#browser()
	mostSensitiveAligner <- calls %>%
		dplyr::select(Id, CallSet, tp) %>%
		group_by(Id, CallSet) %>%
		summarise(tp=sum(tp)) %>%
		ungroup() %>%
		left_join(md, by="Id")

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
	#browser()
	rocby <- function(rocdf, ...) {
		groupingCols <- quos(...)
		rocdf %>%
				dplyr::select(Id, CallSet, !!!groupingCols, QUAL, tp, fp, fn) %>%
				rbind(rocdf %>% dplyr::select(Id, CallSet, !!!groupingCols) %>% distinct(Id, CallSet, !!!groupingCols) %>% mutate(QUAL=max(rocdf$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
				#filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
				group_by(Id, CallSet, !!!groupingCols) %>%
				arrange(desc(QUAL)) %>%
				mutate(events=sum(tp) + sum(fn)) %>%
				mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
				# each QUAL score is a point on the ROC plott
				group_by(Id, CallSet, !!!groupingCols, QUAL, events) %>%
				summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
				# QUAL scores with the same number of tp calls can be merged on the ROC plot
				group_by(Id, CallSet, !!!groupingCols, tp, events) %>%
				summarise(fp=max(fp), fn=max(fn), QUAL=min(QUAL)) %>%
				# subsample along tp and tp+fp axis
				group_by(Id, CallSet, !!!groupingCols) %>%
				dplyr::slice(unique(c(
					1,
					findInterval(seq(0, max(tp), max(tp)/rocSlicePoints), tp),
					findInterval(seq(0, max(tp + fp), max(tp + fp)/rocSlicePoints), tp + fp),
					n()
				))) %>%
				ungroup() %>%
				mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
				left_join(md, by="Id")
	}

	roc <- rocby(calls)
	rocbyrepeat <- rocby(calls, repeatClass)
	# calls$ihomlenBin <- cut(abs(calls$ihomlen), breaks=c(-1000000000, 1, 2, 3, 4, 5, 10, 20, 50, 100, 1000000000), right=FALSE, labels=c("0","1","2", "3", "4", "5-9", "10-19", "20-49", "50-99", "100+"))
	# rocbyrepeat <- rocby(calls, ihomlenBin)
	nominalCalls$snp50bin <- cut(nominalCalls$snp50bp, breaks=c(-1000000000, 1, 2, 3, 6, 10, 50, 100, 1000000000), right=FALSE, labels=c("0","1","2", "3-5", "6-9", "10-49", "50-99", "100+"))
	rocbysnp50 <- rocby(nominalCalls, snp50bin)

	#browser()
	bpErrorDistribution <- bind_rows(calls, nominalCalls) %>%
		filter(tp) %>%
		dplyr::select(Id, CallSet, nominalPosition, bperror) %>%
		group_by(Id, CallSet, nominalPosition, bperror) %>%
		summarize(n=n())
	bpErrorDistribution <- bpErrorDistribution %>%
		left_join(bpErrorDistribution %>% group_by(Id, CallSet, nominalPosition) %>% summarize(count=sum(n))) %>%
		ungroup() %>%
		mutate(rate=n/count) %>%
		dplyr::select(-count) %>%
		left_join(md, by="Id")

	eventSize <- calls %>%
		mutate(Classification = ifelse(tp, "True Positive", ifelse(fp, "False Positive", "False Negative"))) %>%
		dplyr::select(Id, CallSet, svLen, Classification) %>%
		left_join(md, by="Id")

	callsByCaller <- calls %>%
		filter(!fp) %>%
		dplyr::select(-HOMLEN, -ihomlen, -sizeerror, -bperror, -duptp, -fp, -fn, -includeFiltered)

	return(list(mostSensitiveAligner=mostSensitiveAligner,
							callsByEventSize=callsByEventSize,
							roc=roc,
							rocbyrepeat=rocbyrepeat,
							#rocbyihomlen=rocbyihomlen,
							rocbysnp50=rocbysnp50,
							bpErrorDistribution=bpErrorDistribution,
							eventSize=eventSize,
							callsByCaller=callsByCaller))
}
.LoadMinCallsForId <- function(
		ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, eventtypes,
		datadir, metadata, id,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {
	calls <- .CachedLoadCallsForId(datadir, metadata, id,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition)
	if (is.null(calls)) {
		return(NULL)
	}
	calls <- calls %>% mutate(nominalPosition=nominalPosition)
	if (ignore.duplicates) {
		calls <- calls %>% filter(!duptp)
	}
	if (!is.null(eventtypes)) {
		calls <- calls %>% filter(simpleEvent %in% eventtypes)
	}
	if (is.null(metadata$CX_MULTIMAPPING_LOCATIONS)) {
		metadata$CX_MULTIMAPPING_LOCATIONS <- NA_integer_
	}
	calls <- calls %>% mutate(Classification = ifelse(tp, "True Positive", ifelse(fp, "False Positive", "False Negative")))
	return(calls)
}

.CachedLoadCallsForId <- function(datadir, metadata, id,
		ignore.interchromosomal, mineventsize, maxeventsize,
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr, truthgrName,
		grtransform, grtransformName, nominalPosition) {
	if (is.null(truthgrName) || is.na(truthgrName) || truthgrName == "") {
		truthgrName <- GetId((metadata %>% filter(Id==id))$CX_REFERENCE_VCF)
	}
	cachekey <- list(
		datadir=datadir, id=id,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgrName=truthgrName,
		grtransformName=grtransformName, nominalPosition=nominalPosition)
	cachedir <- ".Rcache/Calls"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadCallsForId %s (%s)", id, getChecksum(cachekey)), stderr())
		result <- .LoadCallsForId(
			datadir=datadir, metadata=metadata, id=id,
			ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
			maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, requiredHits=requiredHits, truthgr=truthgr,
			grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.LoadCallsForId <- function(datadir, metadata, id,
		ignore.interchromosomal, mineventsize, maxeventsize,
		# ScoreVariantsFromTruth
		maxgap, sizemargin, ignore.strand, requiredHits, truthgr,
		# .CachedTransformVcf
		grtransform, grtransformName, nominalPosition) {
	if (is.null(truthgr)) {
		truthid <- GetId((metadata %>% filter(Id==id))$CX_REFERENCE_VCF)
		if (is.na(truthid) | is.null(truthid)) {
			stop("Missing truth for ", id)
		}
		truthgr <- .FilteredTransformVcf(datadir=datadir, metadata=metadata, id=truthid, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
	}
	if (is.null(truthgr)) {
		stop("Missing truth for ", id)
	}
	callgr <- .FilteredTransformVcf(datadir=datadir, metadata=metadata, id=id, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
	if (is.null(callgr)) {
		return(NULL)
	}
	mcalls <- NULL
	for (includeFiltered in c(TRUE, FALSE)) {
		calls <- .ScoreVariantsFromTruthVCF(callgr=callgr, truthgr=truthgr, includeFiltered=includeFiltered, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, id=id, requiredHits=requiredHits)
		# can't do straight-forward calls$calls$snp50bp <- callgr$snp50bp if includeFiltered is FALSE
		calls$calls$snp50bp <- callgr$snp50bp[calls$calls$breakendId]
		calls$truth$snp50bp <- truthgr$snp50bp[calls$truth$breakendId]
		mergedcalls <- calls$calls
		if (!is.null(calls$truth)) {
			calls$truth <- calls$truth %>% mutate(duptp=FALSE)
			calls$calls <- calls$calls %>% filter(!tp) # take the truth vcf version of tp calls since it has the actual event size
			mergedcalls <- rbind(calls$calls, calls$truth)
		}
		mcalls <- rbind(mcalls, mergedcalls %>% mutate(CallSet = ifelse(includeFiltered, ALL_CALLS, PASS_CALLS)))
	}
	return(mcalls %>% dplyr::select(-ignore.strand, -maxgap))
}

.cacheableEventTypes <- function(eventtypes) {
	if (!is.null(eventtypes)) {
		# order does not matter; strip names that were confusing the cache
		eventtypes <- as.character(sort(eventtypes))
	}

	return(eventtypes)
}

.CachedLoadCallMatrixForIds <- function(datadir, metadata, ids, eventtypes,
																				ignore.interchromosomal, mineventsize, maxeventsize,
																				maxgap, sizemargin, ignore.strand,
																				# .CachedTransformVcf
																				grtransform, grtransformName, nominalPosition) {
	eventtypes <- .cacheableEventTypes(eventtypes)
	cachekey <- list(datadir, metadata, ids, eventtypes,
							 ignore.interchromosomal, mineventsize, maxeventsize,
							 maxgap, sizemargin, ignore.strand,
							 # .CachedTransformVcf
							 grtransformName, nominalPosition)
	cachedir <- ".Rcache/callmatrix"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".CachedLoadCallMatrixForIds %s %s (%s)", datadir, paste(eventtypes, collapse = ","), getChecksum(cachekey)), stderr())
		result <- .LoadCallMatrixForIds(datadir, metadata, ids, eventtypes,
																		ignore.interchromosomal, mineventsize, maxeventsize,
																		maxgap, sizemargin, ignore.strand,
																		grtransform, grtransformName, nominalPosition)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}

.LoadCallMatrixForIds <- function(datadir, metadata, ids, eventtypes,
		ignore.interchromosomal, mineventsize, maxeventsize,
		maxgap, sizemargin, ignore.strand,
		# .CachedTransformVcf
		grtransform, grtransformName, nominalPosition) {
	eventtypes <- .cacheableEventTypes(eventtypes)
	# include truth in table
	truthids <- GetId((metadata %>% filter(Id %in% ids))$CX_REFERENCE_VCF)
	ids <- unique(c(ids, truthids))
	grlist <- lapply(ids, function(id) {
		gr <- .FilteredTransformVcf(datadir=datadir, metadata=metadata, id=id,
			ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
			grtransform=grtransform, grtransformName=grtransformName, nominalPosition=nominalPosition)
		gr$Id <- id
		names(gr) <- paste0("Id", id, names(gr))
		gr$partner <- paste0("Id", id, gr$partner)
		# Force filters to match on both sides
		gr$FILTER <- ifelse(gr$FILTER %in% c(".", "PASS"), partner(gr)$FILTER, gr$FILTER)

		# Filter on event type
		gr$simpleEvent <- simpleEventType(gr)
		if (!is.null(eventtypes)) {
			gr <- gr[gr$simpleEvent %in% eventtypes]
		}
		return(gr)
	})
	names(grlist) <- ids
	allgr <- GRangesList(grlist)
	allgr <- unlist(allgr, recursive=TRUE, use.names=FALSE)
	allgr$ignore.filtered <- rep(FALSE, length(allgr)) # Crashing here on empty grs for BND

	allgrfiltered <- GRangesList(grlist)
	allgrfiltered <- unlist(allgrfiltered, recursive=TRUE, use.names=FALSE)
	allgrfiltered$ignore.filtered <- rep(TRUE, length(allgrfiltered))
	allgrfiltered <- allgrfiltered[allgrfiltered$FILTER %in% c(".", "PASS")]
	names(allgrfiltered) <- paste0("f", names(allgrfiltered))[seq_along(allgrfiltered)]
	allgrfiltered$partner <- paste0("f", allgrfiltered$partner)[seq_along(allgrfiltered)]

	allgr <- c(allgr, allgrfiltered)
	allgr$CallSet = ifelse(allgr$ignore.filtered, rep(PASS_CALLS, length(allgr)), rep(ALL_CALLS, length(allgr)))

	# initialise to -1 to indicate no match
	for (id in ids) {
		colname <- paste0("Id", id)
		mcols(allgr)[[colname]] <- rep(-1, length(allgr))
		colname <- paste0("fId", id)
		mcols(allgr)[[colname]] <- rep(-1, length(allgr))
	}
	for (ignore.filtered.q in c(TRUE, FALSE)) {
		for (ignore.filtered.s in c(TRUE, FALSE)) {
			for (qid in ids) {
				for (sid in ids) {
					colname <- paste0(ifelse(ignore.filtered.s, "", "f"), "Id", sid)
					qgr <- grlist[[qid]]
					if (ignore.filtered.q) {
						qgr <- qgr[qgr$FILTER %in% c(".", "PASS")]
					}
					sgr <- grlist[[sid]]
					if (ignore.filtered.s) {
						sgr <- sgr[sgr$FILTER %in% c(".", "PASS")]
					}
					key <- list(datadir = datadir, queryid = qid, subjectid = sid, grtransformName = grtransformName, eventtypes = eventtypes)
					if (qid <= sid) {
						cmq <- .CacheMatchingQuals(qgr, sgr, maxgap, sizemargin, ignore.strand, key)
						mcols(allgr[allgr$Id==qid & allgr$ignore.filtered==ignore.filtered.q])[[colname]] <- cmq$bestSubjectQUALforQuery
					} else {
						cmq <- .CacheMatchingQuals(sgr, qgr, maxgap, sizemargin, ignore.strand, key)
						mcols(allgr[allgr$Id==qid & allgr$ignore.filtered==ignore.filtered.q])[[colname]] <- cmq$bestQueryQUALforSubject
					}
				}
			}
		}
	}
	return(allgr)
}
colname_to_CallSet <- function(colname) {
	ifelse(str_detect(colname, "^fId"), ALL_CALLS, PASS_CALLS)
}
colname_to_Id <- function(colname) {
	str_replace(colname, "f?Id", "")
}
IdCallSet_to_colname <- function(id, callset) {
	paste0(ifelse(callset == ALL_CALLS, "fId", "Id"), id)
}
.CacheMatchingQuals <- function(querygr, subjectgr, maxgap, sizemargin, ignore.strand, key) {
	cachekey <- list(maxgap, sizemargin, ignore.strand, key=key,
		# effectively the same as adding ignore.filtered but reuses the cache entry
		queryLength=length(querygr), subjectLength=length(subjectgr))
	cachedir <- ".Rcache/pairwiseQUAL"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".findMatchingQuals %s %s (%s)", key$queryid, key$subjectid, getChecksum(cachekey)), stderr())
		result <- .findMatchingQuals(querygr, subjectgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.findMatchingQuals <- function(querygr, subjectgr, maxgap, sizemargin, ignore.strand, missingQUAL=-1, duplicateHitQUAL=-2) {
		hits <- findBreakpointOverlaps(querygr, subjectgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin)
		hits$queryQUAL <- querygr$QUAL[hits$queryHits]
		hits$subjectQUAL <- subjectgr$QUAL[hits$subjectHits]
		if (any(is.na(hits$queryQUAL))) {
			stop("Missing QUAL")
		}
		if (any(is.na(hits$subjectQUAL))) {
			stop("Missing QUAL")
		}

		bestQueryQUALforSubject <- rep(missingQUAL, length(subjectgr))
		if (is.null(duplicateHitQUAL)) {
			hits <- hits[order(hits$queryQUAL),] # sort by qual so the highest QUAL writes last when doing assignments
			bestQueryQUALforSubject[hits$subjectHits] <- hits$queryQUAL
		} else {
			bestQueryQUALforSubject[hits$subjectHits] <- duplicateHitQUAL
			# take only the highest scoring hit
			hits <- hits[order(-hits$subjectQUAL, hits$queryQUAL),]
			bestHits <- hits[!duplicated(hits$queryHits),]
			bestQueryQUALforSubject[bestHits$subjectHits] <- bestHits$queryQUAL
		}

		bestSubjectQUALforQuery <- rep(missingQUAL, length(querygr))
		if (is.null(duplicateHitQUAL)) {
			hits <- hits[order(hits$subjectQUAL),] # sort by qual so the highest QUAL writes last when doing assignments
			bestSubjectQUALforQuery[hits$queryHits] <- hits$subjectQUAL
		} else {
			bestSubjectQUALforQuery[hits$queryHits] <- duplicateHitQUAL
			# take only the highest scoring hit
			hits <- hits[order(-hits$queryQUAL, hits$subjectQUAL),]
			bestHits <- hits[!duplicated(hits$subjectHits),]
			bestSubjectQUALforQuery[bestHits$queryHits] <- bestHits$subjectQUAL
		}

		return(list(
			bestQueryQUALforSubject=bestQueryQUALforSubject,
			bestSubjectQUALforQuery=bestSubjectQUALforQuery))
}
.FilteredTransformVcf <- function(datadir, metadata, id,
																	ignore.interchromosomal, mineventsize, maxeventsize,
																	grtransform, grtransformName, nominalPosition) {
	cachekey <- list(datadir, id, grtransformName, nominalPosition, ignore.interchromosomal, mineventsize, maxeventsize)
	cachedir <- ".Rcache/fvcfgr"
	assert_that(!is.null(id))
	assert_that(!is.na(id))
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		result <- .CachedTransformVcf(datadir, metadata, id, grtransform, grtransformName, nominalPosition)
		if (!is.null(result)) {
			if (ignore.interchromosomal) {
				result <- result[!is.na(result$svLen)]
			}
			if (!is.null(mineventsize)) {
				result <- result[is.na(result$svLen) | pmax(abs(result$svLen), result$insLen) >= mineventsize]
			}
			if (!is.null(maxeventsize)) {
				result <- result[is.na(result$svLen) | pmax(abs(result$svLen), result$insLen) <= mineventsize]
			}
			# sanity adjustments to filters
			result <- result[names(result) %in% result$partner & result$partner %in% names(result)]

			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	if (!is.null(result) & length(result) == 0) {
		return(NULL)
	}
	return(result)
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
		result <- .TransformVcf(datadir=datadir, id=id, caller=caller, grtransform=grtransform, nominalPosition=nominalPosition, metadata=md)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.TransformVcf <- function(datadir, id, caller, grtransform=NULL, vcftransform=NULL, nominalPosition, metadata) {
	#' Loads structural variant GRanges from the VCFs in the given directory
	filename <- list.files(datadir, pattern=paste0("^", id, ".*.vcf$"), full.names=TRUE)
	if (length(filename) == 0) {
		warning(paste("No vcf for id:", id, "skipping"))
		return(NULL)
	}
	if (length(filename) > 1) {
		browser()
		stop(paste0("Multiple VCFs found for ", id))
	}
	vcf <- readVcf(filename, "")
	if (!is.null(vcftransform)) {
		vcf <- vcftransform(vcf)
	}
	# # homology annotation
	# ihomlendt <- .CachedLoadInexactHomologyAnnotation(datadir, id)
	# if (is.null(ihomlendt)) {
	# 	return(NULL)
	# }
	# # remove breakend suffix and get the longest homology length
	# ihomlendt <- ihomlendt %>%
	# 	dplyr::mutate(sourceId=str_replace(breakpointid, "_bp[0-9]+$", "")) %>%
	# 	dplyr::group_by(sourceId) %>%
	# 	dplyr::summarise(ihomlen=max(ihomlen))
	# ihomlen <- ihomlendt$ihomlen
	# names(ihomlen) <- ihomlendt$sourceId
	vcf <- withqual(vcf, caller)
	gr <- breakpointRanges(vcf, nominalPosition, suffix="_bp")
	if (!is.null(grtransform)) {
		gr <- grtransform(gr, metadata)
	}
	# strip any unpaired breakends
	gr <- gr[gr$partner %in% names(gr)]
	if (length(gr) == 0) {
		return(NULL)
	}
	gr$paramRangeID <- NULL
	gr$REF <- NULL
	gr$ALT <- NULL
	#gr$svtype <- NULL
	#gr$svLen <- NULL
	gr$insSeq <- NULL
	#gr$insLen <- NULL

	# gr$ihomlen <- ihomlen[gr$sourceId]
	# annotate nearby SNP/indel counts
	gr$snp50bp <- 0
	if (!is.null(metadata$CX_SNP_TRUTH)) {
		snpgr <- .CachedRawVcfGRanges(datadir, metadata$CX_SNP_TRUTH)
		seqlevelsStyle(snpgr) <- seqlevelsStyle(gr)
		gr$snp50bp <- countOverlaps(gr, snpgr, maxgap=50)
		if (!nominalPosition) {
			# we want to calculate SNVs in the 50bp window around the nominal
			# position so as to ensure that two callers making the same call
			# report the same snp50 count regardless of the confidence interval
			# of the call
			# the nominal position call set and the confidence interval
			# call sets can differ as grtransform can perform filtering
			# The names should still match since we perform our filtering
			# after we generate the breakpoint GRanges
			nomgr <- gr
			nomgr <- breakpointRanges(vcf, TRUE, suffix="_bp")
			if (!is.null(grtransform)) {
				nomgr <- grtransform(nomgr, metadata)
			}
			nomgr <- nomgr[names(nomgr) %in% names(gr)]
			gr[names(nomgr)]$snp50bp <- countOverlaps(nomgr, snpgr, maxgap=50)
		}
	}
	return(gr)
}
.CachedLoadInexactHomologyAnnotation <- function(datadir, id) {
	cachekey <- list(datadir, id)
	cachedir <- ".Rcache/annotateihomlen"
	assert_that(!is.null(id))
	assert_that(!is.na(id))
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".LoadInexactHomologyAnnotation %s (%s)", id, getChecksum(cachekey)), stderr())
		result <- .LoadInexactHomologyAnnotation(datadir=datadir, id=id)
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
.LoadInexactHomologyAnnotation <- function(datadir, id) {
	#' Loads structural variant GRanges from the VCFs in the given directory
	filename <- list.files(datadir, pattern=paste0("^", id, ".*.vcf$"), full.names=TRUE)
	if (length(filename) == 0) {
		warning(paste("No vcf for id:", id, "skipping"))
		return(NULL)
	}
	if (length(filename) > 1) {
		stop(paste0("Multiple VCFs found for ", id))
	}
	filename <- list.files(datadir, pattern=paste0("^", id, ".*.annimphom.bedpe$"), full.names=TRUE)
	if (length(filename) == 0) {
		warning(paste(".annimphom.bedpe not found for id:", id, "skipping"))
		return(NULL)
	}
	dt <- fread(filename, sep="\t", select=c(7:8))
	names(dt) <- c("breakpointid", "ihomlen")
	return(dt)
}
.CachedRawVcfGRanges <- function(datadir, filename, strip=TRUE) {
	cachekey <- list(datadir, filename, strip)
	cachedir <- ".Rcache/rawvcfgr"
	result <- loadCache(key=cachekey, dirs=cachedir)
	if (is.null(result)) {
		write(sprintf(".CachedRawVcfGRanges %s (%s)", filename, getChecksum(cachekey)), stderr())
		file <- filename
		if (!file.exists(file)) {
			file <- list.files(datadir, pattern=paste0(filename, "$"), full.names=TRUE)
		}
		result <- rowRanges(readVcf(file, ""))
		if (strip) {
			# minimal GRanges object
			names(result) <- NULL
			mcols(result) <- NULL
		}
		if (!is.null(result)) {
			saveCache(result, key=cachekey, dirs=cachedir)
		}
	}
	return(result)
}
