source("global.R")
source("sv_benchmark.R")
source("shinyCache2.R")
library(tidyverse)
library(ggplot2)

rocby <- function(callgr, ..., rocSlicePoints=100, truth_id, ignore.duplicates=TRUE) {
	groupingCols <- quos(...)
	callgr$truth <- mcols(callgr)[[paste0("fId",truth_id)]]
	if (ignore.duplicates) {
		# ignore calls that have a -2 in the truth column (don't touch the truth set though!)
		callgr <- callgr[callgr$Id != truth_id & callgr$truth != -2]
	}
	eventCount <- sum(callgr$Id == truth_id)
	as.data.frame(callgr) %>%
		dplyr::select(Id, CallSet, !!!groupingCols, QUAL, truth) %>%
		group_by(Id, CallSet, !!!groupingCols, QUAL) %>%
		summarise(
			fp=sum(truth < 0),
			tp=sum(truth >= 0)) %>%
		group_by(Id, CallSet, !!!groupingCols) %>%
		arrange(desc(QUAL)) %>%
		mutate(
			fp=cumsum(fp),
			tp=cumsum(tp)) %>%
		# each QUAL score is a point on the ROC plott
		group_by(Id, CallSet, !!!groupingCols, QUAL) %>%
		summarise(tp=max(tp), fp=max(fp)) %>%
		# QUAL scores with the same number of tp calls can be merged on the ROC plot
		group_by(Id, CallSet, !!!groupingCols, tp) %>%
		summarise(fp=max(fp), QUAL=min(QUAL)) %>%
		# subsample along tp and tp+fp axis
		group_by(Id, CallSet, !!!groupingCols) %>%
		dplyr::slice(unique(c(
			1,
			findInterval(seq(0, max(tp), max(tp)/rocSlicePoints), tp),
			findInterval(seq(0, max(tp + fp), max(tp + fp)/rocSlicePoints), tp + fp),
			n()
		))) %>%
		ungroup() %>%
		mutate(
			fn=eventCount-tp,
			precision=tp / (tp + fp),
			fdr=1-precision,
			sens=tp/eventCount)
}

generate_figures <- function(datadir, sample_name, ids, truth_id, truth_name, grtransformName, allow_missing_callers=FALSE) {
	fileprefix <- str_replace(paste(sample_name, truth_name, sep="_"), "[ /]", "_")
	all_ids <- c(truth_id, ids)
	metadata <- LoadCachedMetadata(datadir)
	metadata <- metadata %>% filter(Id %in% all_ids)
	# force truth
	metadata$CX_REFERENCE_VCF <- list.files(datadir, pattern=paste0("^", truth_id, ".*.vcf"))

	missing_callers <- fulldatacallers[!(fulldatacallers %in% StripCallerVersion(metadata$CX_CALLER))]
	if (!allow_missing_callers && length(missing_callers) > 0) {
		stop(paste("Missing metadata for ", missing_callers))
	}

	callgr <- .LoadCallMatrixForIds(
		datadir=datadir,
		metadata=metadata,
		ids=all_ids,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand,
		grtransform=lroptions$grtransform[[grtransformName]],
		grtransformName=grtransformName,
		nominalPosition=nominalPosition)

	callgr$caller_hits_ex_truth <- rowSums(as.matrix(as.data.frame(
		mcols(callgr)[,str_detect(names(mcols(callgr)), "^fId[a-f0-9]+") & !(names(mcols(callgr)) %in% c(paste0("fId", truth_id)))])) != -1)

	missing_callers <- StripCallerVersion((metadata %>% filter(Id %in% unique(callgr$Id) & Id != truth_id & !(CX_CALLER %in% fulldatacallers)))$CX_CALLER)
	if (!allow_missing_callers && length(missing_callers) > 0) {
		stop(paste("Missing VCF records for ", missing_callers))
	}

	caller_colour_scheme <- scale_colour_brewer(palette="Set1")
	prec_recall_ggcommon <-
		# TODO: should we use recall or #tp as x axis scale?
		aes(y=precision, x=tp / 2, colour=Id, linetype=CallSet) +
		caller_colour_scheme

	# Figure 1: prec_recall
	ggplot(rocby(callgr, truth_id=truth_id) %>% filter(Id != truth_id)) +

		geom_line() +
		prec_recall_ggcommon +
		labs(title=paste("Precision-Recall \n", sample_name, truth_name))

	saveplot(paste0(fileprefix, "prec_recall_overall")) #TODO: plot dimensions


	# Figure 2: simulation
	# done by precache.R


	# Figure 3:
	# bar chart
	# precision by #callers
	# MDS/PCA?

	# Figure 4:
	ggplot(rocby(callgr,  truth_id=truth_id) %>% filter(Id != truth_id)) +
		aes(y=precision, x=tp / 2, colour=Id, linetype=CallSet) +
		geom_line() +
		prec_recall_ggcommon +
		labs(title=paste("Precision-Recall by event size\n", sample_name, truth_name))
	saveplot(paste0(fileprefix, "prec_recall_by_event_size"))

}
