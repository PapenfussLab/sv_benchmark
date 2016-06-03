source("sv_benchmark.R")
library(shiny)
library(ggplot2)
library(dplyr)
#' Loads the given call set of simulated data
LoadSimulatedCallSets <- function(datadir, maxgap, ignore.strand, ignore.duplicates, sizemargin) {
	return(LoadCallSets(datadir, maxgap, ignore.strand, ignore.duplicates, sizemargin, vcftransform=function(id, metadata, vcfs) {
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
		}))
}
#' loads the given call set
LoadCallSets <- function(datadir, maxgap, ignore.strand, ignore.duplicates, sizemargin, vcftransform, truthgr=NULL) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	truthhash <- NULL
	if (!is.null(truthgr)) {
	  truthhash <- getChecksum(truthgr)
	}
	key <- list(datadir, maxgap, ignore.strand, ignore.duplicates, sizemargin, truthhash)
	result <- loadCache(key=key, dirs="LoadCallSets")
	if (!is.null(result)) {
		result <- .LoadCallSets(datadir, maxgap, ignore.strand, ignore.duplicates, sizemargin, vcftransform, truthhash)
		saveCache(result, key=key, dirs="LoadCallSets")
	}
	result$key <- key
	setCacheRootPath(cacheroot)
	return(result)
}
.LoadCallSets <- function(datadir, maxgap, ignore.strand, , sizemargin, vcftransform, truthgr) {
	metadata <- LoadMetadata(datadir)
	vcfs <- LoadMinimalSVFromVCF(datadir), metadata=metadata)
	vcfs <- sapply(names(vcfs), vcftransform, metadata, vcfs, simplify=FALSE, USE.NAMES=TRUE)
	mcalls <- NULL
	for (includeFiltered in c(TRUE, FALSE)) {
		calls <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=includeFiltered, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand, truthgr=truthgr)
		calls$truth <- calls$truth %>% mutate(duptp=FALSE)
		calls$calls <- calls$calls%>% filter(!tp) # take the truth vcf version of tp calls since it has the actual event size
		mergedcalls <- rbind(calls$calls, calls$truth)
		mergedcalls$CallSet <- ifelse(includeFiltered ? "High & Low confidence" : "High confidence only")
		mcalls <- rbind(mcalls, mergedcalls)
	}
	return(list(metadata=metadata, calls=mcalls))
}
LoadGraphDataFrames <- function(datadir, callset, ignore.duplicates) {
	cacheroot <- getCacheRootPath()
	setCacheRootPath(datadir)
	key <- list(callset$key, ignore.duplicates)
	result <- loadCache(key=key, dirs="LoadGraphDataFrames")
	if (!is.null(result)) {
		result <- .LoadGraphDataFrames(callset, ignore.duplicates)
		saveCache(result, key=key, dirs="LoadGraphDataFrames")
	}
	setCacheRootPath(cacheroot)
	return(result)
}
.LoadGraphDataFrames <- function(callset, ignore.duplicates, mineventsize=NULL, maxeventsize=NULL) {
	if (ignore.duplicates) {
		callset$calls <- callset$calls %>% filter(!duptp)
	}
	# TODO event size filters
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

shinyServer(function(input, output) {
  output$distPlot <- renderPlot({

    # draw the histogram with the specified number of bins
    hist(x, breaks = 1, col = 'darkgray', border = 'white')
  })
})
