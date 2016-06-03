library(GenomicRanges)
library(StructuralVariantAnnotation) #install_github("d-cameron/StructuralVariantAnnotation")
library(testthat)
library(stringr)
library(dplyr)
library(tidyr)
library(R.cache)

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/projects/sv_benchmark/", "~/projects/sv_benchmark/")

GetId <- function(filenames) {
	cf <- as.character(filenames)
	if (length(cf) == 0) {
		return(character(0))
	} else if (all(is.na(cf))) {
		return(cf)
	} else {
		return(str_match(basename(as.character(filenames)), "([^.]+)(\\..*)*$")[,2])
	}
}
# Load metadata into a dataframe
LoadMetadata <- function(directory) {
	write("Loading metadata", stderr())
	filenames <- list.files(directory, pattern="*.metadata$", full.names=TRUE)
	zeroSizeFiles = file.info(filenames)$size == 0
	if (any(zeroSizeFiles)) {
		warning(paste("Skipping files", filenames[zeroSizeFiles], "as they have 0 size."))
		filenames <- filenames[!zeroSizeFiles]
	}
	#metadata <- foreach (filename=filenames, .export=c("GetId"), .combine=rbind) %dopar% {
	metadata <- lapply(filenames, function(filename) {
		md <- read.csv(filename, header=FALSE, sep="=", quote = "\"'", col.names=c("CX", "V"))
		md$File <- filename
		md$Id <- GetId(filename)
		md
	})
	metadata <- do.call(rbind, metadata)
	metadata <- data.frame(lapply(metadata, as.character), stringsAsFactors=FALSE)
	metadata <- spread(metadata, CX, V)
	# convert data from older format
	if (!is.null(metadata$CX_ALIGNER_SOFTCLIP)) {
		metadata$CX_ALIGNER_MODE <- metadata$CX_ALIGNER_MODE %na% ifelse(metadata$CX_ALIGNER_SOFTCLIP == 1, "local", "global")
	}
	# transform known numeric data to expected type
	if (!is.null(metadata$CX_READ_FRAGMENT_LENGTH)) {
		metadata$CX_READ_FRAGMENT_LENGTH <- as.numeric(as.character(metadata$CX_READ_FRAGMENT_LENGTH))
	}
	if (!is.null(metadata$CX_READ_LENGTH)) {
		metadata$CX_READ_LENGTH <- as.numeric(as.character(metadata$CX_READ_LENGTH))
	}
	if (!is.null(metadata$CX_READ_DEPTH)) {
		metadata$CX_READ_DEPTH <- as.numeric(as.character(metadata$CX_READ_DEPTH))
	}
	rownames(metadata) <- metadata$Id
	write(paste(nrow(metadata), "metadata files loaded"), stderr())
	return(metadata)
}
# infer proxy quality scores for ROC purposes based on strength of support
withqual <- function(vcf, caller) {
	if (is.null(rowRanges(vcf)$QUAL)) {
		rowRanges(vcf)$QUAL <- NA_real_
	}
	if (!is.na(caller) && !is.null(caller) && all(is.na(rowRanges(vcf)$QUAL))) {
		caller <- str_extract(caller, "^[^/]+") # strip version
		# use total read support as a qual proxy
		if (caller %in% c("delly")) {
			rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$PE), 0, info(vcf)$PE) + ifelse(is.na(info(vcf)$SR), 0, info(vcf)$SR)
		} else if (caller %in% c("crest")) {
			rowRanges(vcf)$QUAL <- ifelse(is.na(info(vcf)$right_softclipped_read_count), 0, info(vcf)$right_softclipped_read_count) + ifelse(is.na(info(vcf)$left_softclipped_read_count), 0, info(vcf)$left_softclipped_read_count)
		} else if (caller %in% c("pindel")) {
			rowRanges(vcf)$QUAL <- geno(vcf)$AD[,1,2]
		} else if (caller %in% c("lumpy")) {
			rowRanges(vcf)$QUAL <- unlist(info(vcf)$SU)
		} else if (caller %in% c("cortex")) {
			rowRanges(vcf)$QUAL <- geno(vcf)$COV[,1,1]
		}
	}
	if (any(is.na(rowRanges(vcf)$QUAL))) {
		#if (is.null(caller) && is.na(caller)) {
		warning(paste("Missing QUAL scores for", caller))
	}
	return(vcf)
}

StripCallerVersion <- function(caller) {
	caller <- paste0(str_extract(caller, "^([^/]+)"), str_match(caller, "^([^/]+)\\/[^/]+(/[^/]+)?")[,3] %na% "") %na% caller
	if (any(caller=="gridss")) {
	  caller <- relevel(factor(caller), "gridss")
	}
	return(caller)
}
PrettyVariants <- function(x) {
	x[x=="hetDEL"] <- "Deletion"
	x[x=="hetINS"] <- "Insertion"
	x[x=="hetDUP"] <- "Tandem Duplication"
	x[x=="hetINV"] <- "Inversion"
	x[x=="hetBP"] <- "Breakpoint"
	x[x=="hetBP_SINE"] <- "SINE/ALU Breakpoint"
	return(x)
}
#' Loads a minimal structural variant GRanges from the VCF
LoadMinimalSVs <- function(filename, caller, transform=NULL) {
  key <- list(filename, caller, transform)
  gr <- loadCache(key=key, dirs="LoadMinimalSVs")
  if (is.null(gr)) {
  	gr <- .LoadMinimalSVs(filename, caller, transform)
  	saveCache(gr, key=key, dirs="LoadMinimalSVs")
  }
	return(gr)
}
.LoadMinimalSVs <- function(filename, caller, transform=NULL) {
  vcf <- readVcf(filename, "hg19")
  if (!is.null(transform)) {
    vcf <- transform(vcf)
  }
  vcf <- withqual(vcf, caller)
  gr <- breakpointRanges(vcf)
  gr$paramRangeID <- NULL
  gr$REF <- NULL
  gr$ALT <- NULL
  #gr$svtype <- NULL
  #gr$svLen <- NULL
  gr$insSeq <- NULL
  #gr$insLen <- NULL
  return(gr)
}
#vcf <- readVcf("C:/dev/sv_benchmark/data.aligner/5afa7ffdf2cc32602476526d5b477c5c.vcf", "hg19")
#' Loads structural variant GRanges from the VCFs in the given directory
LoadMinimalSVFromVCF <- function(directory, pattern="*.vcf$", metadata=NULL, existingList=NULL, transform=NULL) {
	write("Loading VCFs", stderr())
	filenames <- list.files(directory, pattern=pattern, full.names=TRUE)
	zeroSizeFiles = file.info(filenames)$size == 0
	if (any(zeroSizeFiles)) {
		write(paste("Skipping file", filenames[zeroSizeFiles], "due to 0 size.\n"))
		warning(paste("Skipping files", paste(filenames[zeroSizeFiles]), "due to 0 size.\n"))
		filenames <- filenames[!zeroSizeFiles]
	}
	# exclude already loaded VCFs
	filenames <- filenames[!(GetId(filenames) %in% names(existingList))]
	# exclude VCFs without metadata
	filenames <- filenames[GetId(filenames) %in% metadata$Id]
	# only load VCFS that have metadata
	#if (!is.null(metadata)) {
	#	filenames <- filenames[GetId(filenames) %in% metadata$Id]
	#}
	#vcfs <- foreach (filename=filenames, .packages="VariantAnnotation") %dopar% { # Parallel load of VCFs
	grlist <- lapply(filenames, function(filename) {
		write(paste0("Loading ", filename), stderr())
		caller <- NULL
		if (!is.null(metadata)) {
			caller <- metadata$CX_CALLER[metadata$Id == GetId(filename)]
		}
		gr <- LoadMinimalSVs(filename, caller, transform)
		gr$Id <- rep(GetId(filename), length(gr))
		return (gr)
	})
	names(grlist) <- GetId(filenames)
	grlist[sapply(grlist, is.null)] <- NULL # Remove NULL VCFs list
	write(paste("Loaded", length(grlist), "VCFs"), stderr())
	return(c(existingList, grlist))
}
.interval_distance <- function(s1, e1, s2, e2) {
	return (ifelse(s2 >= s1 & s2 <= e1, 0,
					ifelse(s1 >= s2 & s1 <= e2, 0,
					ifelse(s1 < s2, s2 - e1, s1 - e2))))
}
.distance <- function(r1, r2) {
	return(data.frame(
		min=pmax(0, pmax(start(r1), start(r2)) - pmin(end(r1), end(r2))),
		max=pmax(end(r2) - start(r1), end(r1) - start(r2))))
}
#' Finds matchings breakpoints
#'
#' @param sizemargin error margin in allowable size
#' @param restrictMarginToSizeMultiple size restriction multiplier on event size.
#' The default value of 0.5 requires that the breakpoint positions can be off by
#' at maximum, half the event size. This ensures that small deletion do actually
#' overlap at least one base pair.
#'
#'@export
findMatchingBreakpoints <- function(query, subject, maxgap=0L, ignore.strand=FALSE, sizemargin=0.25, restrictMarginToSizeMultiple=0.5) {
  hits <- findBreakpointOverlaps(query, subject, maxgap=maxgap, ignore.strand=ignore.strand)
  # take into account confidence intervals when calculating event size
  callwidth <- .distance(query, partner(query))
  truthwidth <- .distance(subject, partner(subject))
  callsize <- callwidth + (query$insLen %na% 0)
  truthsize <- truthwidth + (subject$insLen %na% 0)
  hits$sizeerror <- .distance(IRanges(start=callsize$min[hits$queryHits], end=callsize$max[hits$queryHits]),
                              IRanges(start=truthsize$min[hits$subjectHits], end=truthsize$max[hits$subjectHits]))$min
  # event sizes must be within sizemargin
  hits <- hits[hits$sizeerror - 2 < sizemargin * pmin(callsize$max[hits$queryHits], truthsize$max[hits$subjectHits]),]
  # further restrict breakpoint positions for small events
  hits$localbperror <- .distance(query[hits$queryHits], subject[hits$subjectHits])$min
  hits$remotebperror <- .distance(partner(query)[hits$queryHits], partner(subject)[hits$subjectHits])$min
  if (!is.null(restrictMarginToSizeMultiple)) {
    allowablePositionError <- (pmin(callsize$max[hits$queryHits], truthsize$max[hits$subjectHits]) * restrictMarginToSizeMultiple + 2)
    hits <- hits[hits$localbperror <= allowablePositionError & hits$remotebperror <= allowablePositionError, ]
  }
  return(hits)
}
#' Finds pairs of query breakpoints that span a matching subject breakpoint
findSpanningBreakpoints <- function(query, subject, maxgap=0L, ignore.strand=FALSE, sizemargin=0.25, restrictMarginToSizeMultiple=0.5, maxSpanningFragmentSize,
                                    matchDirection=TRUE) {
  # must be larger than min fragment size (if not, small indel + large event get matched to large event)
  query <- query[is.na(query$svLen) | abs(query$svLen) >= maxSpanningFragmentSize]
  # <A1-------A2>   <B1-----B2>
  #             |-d-|
  # find A2-B1 matches
  # d must be < maxSpanningFragmentSize
  spanningHits <- as.data.frame(findOverlaps(query, query, maxgap=maxSpanningFragmentSize, ignore.strand=TRUE))
  # removes self-intersections and duplication due to picking up both sides
  spanningHits <- spanningHits[spanningHits$subjectHits < spanningHits$queryHits,]
  a2 <- query[spanningHits$subjectHits]
  b1 <- query[spanningHits$queryHits]
  spanningHits$fragmentSize <- abs((start(a2) + end(a2)) / 2 - (start(b1) + end(b1)) / 2)
  # directions must be appropriate for a small fragment
  if (matchDirection) {
    spanningHits <- spanningHits[ifelse(paste0(strand(a2), strand(b1)) == "+-", start(a2) < end(b1),
                                 ifelse(paste0(strand(b1), strand(a2)) == "+-", start(b1) < end(a2), FALSE)),]
  }
  # generate A1-B2 breakpoint gr
  a1 <- partner(query)[spanningHits$subjectHits]
  b2 <- partner(query)[spanningHits$queryHits]
  spanninggr <- c(a1, b2)
  spanninggr$partner <- c(names(b2), names(a1))
  spanninggr$fragmentSize <- spanningHits$fragmentSize
  # spanning against subject
  hits <- findMatchingBreakpoints(spanninggr, subject, maxgap, ignore.strand, sizemargin, restrictMarginToSizeMultiple)
  hits$localBreakend <- names(spanninggr)[hits$queryHits]
  hits$remoteBreakend <- names(partner(spanninggr))[hits$queryHits]
  hits$fragmentSize <- spanninggr$fragmentSize[hits$queryHits]
  hits$queryHits <- NULL
  return(hits)
}
#' @param truthhash hash of truthgr if it is not the 'natural' truth to compare to
ScoreVariantsFromTruthVCF <- function(callgr, truthgr, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin=0.25, id=NULL, truthhash=NULL) {
  if (length(callgr) == 0) {
    return(.ScoreVariantsFromTruthVCF(callgr, truthgr, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin=0.25, id %null% NA_character_))
  }
  id <- id %null% callgr$Id[1]
  key <- list(includeFiltered, maxgap, ignore.strand, sizemargin, id, truthhash)
  result <- loadCache(key=key, dirs="ScoreVariantsFromTruthVCF")
  if (is.null(result)) {
    result <- .ScoreVariantsFromTruthVCF(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id)
    saveCache(result, key=key, dirs="ScoreVariantsFromTruthVCF")
  }
  return(result)
}
.ScoreVariantsFromTruthVCF <- function(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id) {
	if (!includeFiltered) {
		callgr <- callgr[callgr$FILTER %in% c("PASS", "."),]
	}
	if (is.null(truthgr)) {
		stop(paste("Missing truth ", truthid, " for ", id))
	}
	hits <- findMatchingBreakpoints(callgr, truthgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin)

	hits$QUAL <- callgr$QUAL[hits$queryHits]
	hits <- hits[order(-hits$QUAL),]

	calldf <- as.data.frame(callgr) %>%
		dplyr::select(Id, QUAL, svLen, insLen, vcfId) %>%
		mutate(tp=FALSE, duptp=FALSE, fp=FALSE, fn=FALSE, sizeerror=NA, bperror=NA) %>%
		mutate(
			includeFiltered=includeFiltered,
			maxgap=maxgap,
			ignore.strand=ignore.strand)
	calldf$tp[hits$queryHits] <- TRUE
	calldf$duptp[hits[duplicated(hits$subjectHits),]$queryHits] <- TRUE
	calldf$fp <- !calldf$tp
	truthdf <- as.data.frame(truthgr) %>%
		dplyr::select(svLen, insLen, vcfId) %>%
		mutate(Id=id, QUAL=0, tp=FALSE, fp=FALSE, fn=FALSE, sizeerror=NA, bperror=NA) %>%
		mutate(
			includeFiltered=includeFiltered,
			maxgap=maxgap,
			ignore.strand=ignore.strand)
	truthdf$tp[hits$subjectHits] <- TRUE
	truthdf$fn <- !truthdf$tp

	truthdf$QUAL[hits$subjectHits] <- hits$QUAL
	calldf$bperror[hits$queryHits] <- hits$localbperror
	truthdf$bperror[hits$subjectHits] <- hits$localbperror
	calldf$sizeerror[hits$queryHits] <- hits$sizeerror
	truthdf$sizeerror[hits$subjectHits] <- hits$sizeerror
	return(list(calls=calldf, truth=truthdf))
}

ScoreVariantsFromTruth <- function(vcfs, metadata, includeFiltered=FALSE, maxgap, ignore.strand, sizemargin=0.25, truthgr=NULL) {
	ids <- metadata$Id[!is.na(metadata$CX_CALLER) & metadata$Id %in% names(vcfs)]
	truthhash <- NULL
	if (!is.null(truthgr)) {
	  truthhash <- getChecksum(truthgr)
	}
	scores <- lapply(ids, function(id) {
		write(paste0("ScoreVariantsFromTruth ", id), stderr())
		callgr <- vcfs[[id]]
		if (is.null(truthgr)) {
			truthid <- GetId((metadata %>% filter(Id==id))$CX_REFERENCE_VCF)
			truthgr <- vcfs[[truthid]]
		}
		if (is.null(truthgr)) {
			stop("Missing truth for ", id)
		}
		return(ScoreVariantsFromTruthVCF(callgr, truthgr, includeFiltered, maxgap, ignore.strand, sizemargin, id, truthhash=truthhash))
	})
	return(list(
		calls=rbind_all(lapply(scores, function(x) x$calls)),
		truth=rbind_all(lapply(scores, function(x) x$truth))))
}

#' subsets the breakpoints to only include breakpoints in which both breakends
#' occur within the specified bed regions
subsetbed <- function(gr, bed, maxgap) {
	gr <- subsetByOverlaps(gr, bed, maxgap=maxgap, ignore.strand=TRUE)
	gr <- gr[gr$partner %in% names(gr)]
	return(gr)
}

LoadVCFs <- function(datadir, result) {
	key <- list(datadir)
	if (!is.null(result) && getChecksum(result$key) == getChecksum(key)) {
		return(result)
	}
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	result <- loadCache(key=key, dirs="LoadVCFs")
	if (is.null(result)) {
		result <- .LoadVCFs(datadir)
		saveCache(result, key=key, dirs="LoadVCFs")
	}
	result$key <- key
	setCacheRootPath(cacheroot)
	return(result)
}
.LoadVCFs <- function(datadir) {
	metadata <- LoadMetadata(datadir)
	vcfs <- LoadMinimalSVFromVCF(datadir, metadata=metadata)
	return(list(metadata=metadata, vcfs=vcfs, datadir=datadir))
}
#' Loads the given call set of simulated data
LoadSimulatedCallSets <- function(mvcfs, maxgap, ignore.strand, sizemargin, result) {
	return(LoadCallSets(mvcfs, maxgap, ignore.strand, sizemargin, vcftransform=function(id, metadata, vcfs) {
			gr <- vcfs[[id]]
			if (!is.na((metadata %>% filter(Id == id))$CX_CALLER)) {
				# only looking at intrachromosomal calls
				gr <- gr[!is.na(gr$svLen),]
				if (str_detect((metadata %>% filter(Id == id))$CX_REFERENCE_VCF_VARIANTS, "BP")) {
					# filtering small events calls to remove spurious indels caused by sequence homology around breakpoints
					gr <- gr[abs(gr$svLen) >= 50,]
				}
			}
			return(gr)
		}, result=result))
}
#' loads the given call set
LoadCallSets <- function(mvcfs, maxgap, ignore.strand, sizemargin, vcftransform, truthgr=NULL, result=NULL) {
	truthhash <- NULL
	if (!is.null(truthgr)) {
	  truthhash <- getChecksum(truthgr)
	}
	key <- list(mvcfs$datadir, maxgap, ignore.strand, sizemargin, truthhash)
	if (!is.null(result) && getChecksum(result$key) == getChecksum(key)) {
		return(result)
	}
	cacheroot <- getCacheRootPath()
	setCacheRootPath(mvcfs$datadir)
	result <- loadCache(key=key, dirs="LoadCallSets")
	if (is.null(result)) {
		result <- .LoadCallSets(mvcfs, maxgap, ignore.strand, sizemargin, vcftransform, truthhash)
		saveCache(result, key=key, dirs="LoadCallSets")
	}
	result$key <- key
	setCacheRootPath(cacheroot)
	return(result)
}
.LoadCallSets <- function(mvcfs, maxgap, ignore.strand, sizemargin, vcftransform, truthgr) {
	metadata <- mvcfs$metadata
	vcfs <- mvcfs$vcfs
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
	return(list(metadata=metadata, calls=mcalls, datadir=mvcfs$datadir))
}
LoadGraphDataFrames <- function(callset, ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize, result=NULL) {
	key <- list(callset$key, ignore.duplicates, ignore.interchromosomal, mineventsize, maxeventsize)
	if (!is.null(result) && getChecksum(result$key) == getChecksum(key)) {
		return(result)
	}
	cacheroot <- getCacheRootPath()
	setCacheRootPath(callset$datadir)
	result <- loadCache(key=key, dirs="LoadGraphDataFrames")
	if (is.null(result)) {
		result <- .LoadGraphDataFrames(callset, ignore.duplicates)
		saveCache(result, key=key, dirs="LoadGraphDataFrames")
	}
	result$key <- key
	setCacheRootPath(cacheroot)
	return(result)
}
.LoadGraphDataFrames <- function(callset, ignore.duplicates, ignore.interchromosomal=TRUE, mineventsize=NULL, maxeventsize=NULL) {
	if (ignore.duplicates) {
		callset$calls <- callset$calls %>% filter(!duptp)
	}
	if (ignore.interchromosomal) {
		callset$calls <- callset$calls %>% filter(!is.na(svLen))
	}
	if (!is.null(mineventsize)) {
		callset$calls <- callset$calls %>% filter(svLen >= mineventsize)
	}
	if (!is.null(maxeventsize)) {
		callset$calls <- callset$calls %>% filter(svLen <= maxeventsize)
	}
	mostSensitiveAligner <- callset$calls %>%
		select(Id, CallSet, tp) %>%
		group_by(Id, CallSet) %>%
		summarise(tp=sum(tp)) %>%
		ungroup() %>%
		arrange(desc(tp)) %>%
		left_join(callset$metadata) %>%
		distinct(CallSet, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF) %>%
		select(Id, CallSet)
	callsByEventSize <- callset$calls %>%
	  #filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
	  filter(!fp) %>%
	  filter(Id %in% callset$metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	  group_by(Id, CallSet, svLen) %>%
	  summarise(sens=sum(tp)/sum(tp+fn)) %>%
	  ungroup() %>%
	  left_join(callset$metadata) %>%
	  mutate(eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS))
	roc <- callset$calls %>%
	  select(Id, CallSet, QUAL, tp, fp, fn) %>%
	  rbind(mcalls %>% select(Id, CallSet) %>% distinct(Id, CallSet) %>% mutate(QUAL=max(mcalls$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
	  #filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
	  group_by(Id, CallSet) %>%
	  arrange(desc(QUAL)) %>%
	  mutate(events=sum(tp) + sum(fn)) %>%
	  mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
	  group_by(Id, CallSet, QUAL, events) %>%
	  summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
	  ungroup() %>%
	  mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
	  left_join(callset$metadata) %>%
	  mutate(eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS))
	return(list(mostSensitiveAligner=mostSensitiveAligner, callsByEventSize=callsByEventSize, roc=roc))
}
