source("global.R")
source("sv_benchmark.R")
source("shinyCache2.R")
library(tidyverse)
library(ggplot2)

rocby <- function(callgr, ..., rocSlicePoints=100, truth_id, ignore.duplicates=TRUE) {
	groupingCols <- quos(...)
	callgr$truthQUAL <- mcols(callgr)[[paste0("fId",truth_id)]]
	if (!ignore.duplicates) {
		rocdf <- callgr %>%
			as.data.frame() %>%
			dplyr::select(Id, CallSet, !!!groupingCols, QUAL, truthQUAL) %>%
			mutate(fp=truthQUAL < 0, tp=truthQUAL > 0)
	} else {
		# ignore calls that have a -2 in the truth column (don't touch the truth set though!)
		callgr <- callgr[callgr$Id == truth_id | callgr$truthQUAL != -2]
		truthhitsdf <- callgr %>%
			as.data.frame() %>%
			filter(Id == truth_id & CallSet == "High & Low confidence") %>%
			dplyr::select(!!!groupingCols, dplyr::matches("f?Id.+", ignore.case=FALSE)) %>%
			gather(key="Id_CallSet", value="QUAL", dplyr::matches("f?Id.+", ignore.case=FALSE)) %>%
			filter() %>%
			mutate(
				fp=0,
				tp=QUAL >= 0,
				fn=QUAL < 0,
				CallSet=ifelse(str_detect(Id_CallSet, "^fId"),"High & Low confidence", "High confidence only"),
				Id=str_replace(Id_CallSet, "f?Id", "")) %>%
			dplyr::select(-Id_CallSet)
		rocdf <- callgr %>%
			as.data.frame() %>%
			mutate(fp=1, tp=0, fn=0) %>%
			filter(Id != truth_id & truthQUAL < 0) %>%
			dplyr::select(Id, CallSet, !!!groupingCols, QUAL, fp, tp, fn) %>%
			rbind(truthhitsdf)
	}
	return (rocdf %>%
		#filter(Id != truth_id) %>%
		group_by(Id, CallSet, !!!groupingCols, QUAL) %>%
		summarise(fp=sum(fp), tp=sum(tp)) %>%
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
		#TODO: fn and sens need eventCount using group_by(Id, CallSet, !!!groupingCols)
		mutate(
			precision=tp / (tp + fp),
			fdr=1-precision))
}

generate_figures <- function(datadir, sample_name, ids, truth_id, truth_name, grtransformName, allow_missing_callers=FALSE) {
	fileprefix <- str_replace(paste(sample_name, truth_name, sep="_"), "[ /]", "_")
	all_ids <- c(truth_id, ids)
	metadata <- LoadCachedMetadata(datadir)
	metadata <- metadata %>% filter(Id %in% all_ids)
	# force truth
	metadata$CX_REFERENCE_VCF <- list.files(datadir, pattern=paste0("^", truth_id, ".*.vcf$"))

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


	callgr$longreadhits <-

	callgr$caller_hits_ex_truth <- rowSums(as.matrix(as.data.frame(
		mcols(callgr)[,str_detect(names(mcols(callgr)), "^fId[a-f0-9]+") & !(names(mcols(callgr)) %in% c(paste0("fId", truth_id)))])) != -1)
	callgr$simpleEvent <- simpleEventType(callgr)

	missing_callers <- StripCallerVersion((metadata %>% filter(Id %in% unique(callgr$Id) & Id != truth_id & !(CX_CALLER %in% fulldatacallers)))$CX_CALLER)
	if (!allow_missing_callers && length(missing_callers) > 0) {
		stop(paste("Missing VCF records for ", missing_callers))
	}

	caller_colour_scheme <- scale_colour_brewer(palette="Set1")

	# Figure 1: prec_recall
	ggplot(rocby(callgr, simpleEvent, truth_id=truth_id) %>%
			filter(Id != truth_id) %>%
			left_join(metadata)) +
		# TODO: should we use recall or #tp as x axis scale?
		aes(y=precision, x=tp / 2, colour=StripCallerVersion(CX_CALLER), linetype=CallSet) +
		geom_line() +
		facet_wrap(~ simpleEvent, scale="free") +
		caller_colour_scheme +
		labs(title=paste("Precision-Recall \n", sample_name, truth_name))

	saveplot(paste0(fileprefix, "prec_recall_overall")) #TODO: plot dimensions


	# Figure 2: simulation
	# done by precache.R


	# Figure 3:
	# bar chart
	# precision by #callers
	# MDS/PCA?

	# Figure 4:
	callgr$snp50bpbin <- cut(callgr$snp50bp, breaks=c(0, 1, 2, 3, 4, 5, 1000), labels=c(0,1,2,3,4,"5+"), right=FALSE)
	ggplot(rocby(callgr, snp50bpbin, truth_id=truth_id)) +
		aes(y=precision, x=tp / 2, colour=Id, linetype=CallSet) +
		geom_line() +
		caller_colour_scheme +
		facet_grid(snp50bpbin ~ ., scales="free") +
		labs(title=paste("Precision-Recall by flanking SNV/indels within 50bp of SV call\n", sample_name, truth_name))
	saveplot(paste0(fileprefix, "prec_recall_by_snvindel"))

	callgr$eventSizeBin <- cut(abs(callgr$svLen), breaks=c(0, 100, 200, 300, 500, 1000, 1000000000), labels=c("50-99", "100-199", "200-299", "300-499", "500-999", "1000+"), right=FALSE)
	ggplot(rocby(callgr, simpleEvent, eventSizeBin, truth_id=truth_id) %>%
			#filter(Id != truth_id) %>%
			left_join(metadata)) +
		aes(y=precision, x=tp / 2, colour=StripCallerVersion(CX_CALLER), linetype=CallSet) +
		geom_line() +
		caller_colour_scheme +
		facet_grid(simpleEvent ~ eventSizeBin, scales="free") +
		labs(title=paste("Precision-Recall by event size and type\n", sample_name, truth_name),
			colour="Caller")

	ggplot(rocby(callgr, repeatClass, simpleEvent, truth_id=truth_id) %>%
			filter(Id != truth_id) %>%
			left_join(metadata)) +
		aes(y=precision, x=tp / 2, colour=StripCallerVersion(CX_CALLER), linetype=CallSet) +
		geom_line() +
		caller_colour_scheme +
		facet_wrap(simpleEvent ~ repeatClass, scales="free") +
		scale_y_continuous(limits=c(0,1)) +
		labs(title=paste("Precision-Recall by event size and type\n", sample_name, truth_name),
			colour="Caller")

	callgr$trf <- overlapsAny(callgr, grtrf[[1]], type="any")
	callgr$repeatAnn <- ifelse(callgr$repeatClass == "" & callgr$trf, "TRF", callgr$repeatClass)
	ggplot(rocby(callgr, repeatAnn, truth_id=truth_id) %>% filter(Id != truth_id) %>% left_join(metadata)) +
		aes(y=precision, x=tp / 2, colour=StripCallerVersion(CX_CALLER), linetype=CallSet) +
		geom_line() +
		caller_colour_scheme +
		facet_wrap( ~ repeatAnn, scales="free") +
		scale_y_continuous(limits=c(0,1)) +
		labs(title=paste("Precision-Recall by presence of tandem repeat at breakpoint\n", sample_name, truth_name), colour="Caller")
}




