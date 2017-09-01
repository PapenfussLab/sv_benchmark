source("global.R")
source("sv_benchmark.R")
source("shinyCache2.R")
library(tidyverse)
library(ggplot2)

ignore.interchromosomal <- TRUE
mineventsize <- 51
maxeventsize <- NULL #TODO: what max event size should we use? Need one for chm since truth only has calls under a certain size
maxgap <- 200
sizemargin <- 0.25
ignore.strand <- TRUE
nominalPosition <- FALSE
grtransformName <- "None"

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

generate_figures <- function(datadir, sample_name, ids, truth_id, truth_name, allow_missing_callers=FALSE) {
	fileprefix <- str_replace(paste(sample_name, truth_name, sep="_"), "[ /]", "_")
	all_ids <- c(truth_id, ids)
	metadata <- LoadCachedMetadata(datadir)
	metadata <- metadata %>% filter(Id %in% all_ids)
	if (is.null(metadata$CX_REFERENCE_VCF)) {
		metadata$CX_REFERENCE_VCF <- list.files(datadir, pattern=paste0("^", truth_id, ".*.vcf"))
	}

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

	callgr$caller_hits_ex_truth <- rowSums(as.matrix(as.data.frame(mcols(callgr)[,
	str_detect(names(mcols(callgr)), "^fId[a-f0-9]+") & !(names(mcols(callgr)) %in% c(
		paste0("fId", truth_id)))])) != -1)


	missing_callers <- StripCallerVersion((metadata %>% filter(Id %in% unique(callgr$Id) & Id != truth_id & !(CX_CALLER %in% fulldatacallers)))$CX_CALLER)
	if (!allow_missing_callers && length(missing_callers) > 0) {
		stop(paste("Missing VCF records for ", missing_callers))
	}

	# Figure 1: prec_recall
	ggplot(rocby(callgr, truth_id=truth_id) %>% filter(Id != truth_id)) +
		aes(y=precision, x=tp, colour=Id, linetype=CallSet) + #TODO: rescale tps to breakpoint counts
		geom_line() +
		scale_colour_brewer(palette="Set1") #TODO: colour scheme
		labs(title=paste("Precision-Recall \n", sample_name, truth_name))

	saveplot(paste0(fileprefix, "prec_recall_overall")) #TODO: plot dimensions
}

# NA12878
na12878_truth <- c(
		# TODO: which truth sets are these?
		"0000"="00000000000000000000000000000000",
		"0001"="00000000000000000000000000000001",
		"0002"="00000000000000000000000000000002",
		"Parikh et al"="00000000000000000000000000000003")
for (i in seq_along(na12878_truth)) {
	datadir <- "../data.na12878"
	sample_name <- "NA12878"
	ids <- c(
		#"2f4d15f5f6e428fb91e12dac971571b9", #bcftools
		#"5cdeeb9fd824f642f0ebee04627adf6e", gridss 1.2.1
		"6ae03359fcf0a39aa236a0aaea7ea915",
		"8aaf2886ffe782e666661d6b890b330a",
		"80dd0c2aa34964f330702ee3d0bfa53f",
		#"1139a1e6ef2d598098fe5c2ff609052a", manta/0.29.6
		"26511afa5055a939fdf99a2ac24938cc",
		"a6cbc4fc76337d871ef899d51fbae0d9",
		"b6aeda0d81ff8c839e381b57898dc3d8",
		#"e77293fb821458c3fd04c58ea88b7b16", #gasv
		#"ee44a2b21004af249f115bb5d9508ceb", gridss 1.3.0
		#"f75b305c5579449e347c1b87832e65d8", #gasv
		"f849b9493f2dd5dc20b6b7e49e1c89d7",
		"fa5f9a52d3b67a206cb6e581053a9269")
	truth_id <- na12878_truth[i]
	truth_name <- names(na12878_truth)[i]
	generate_figures(datadir, sample_name, ids, truth_id, truth_name)
}
#chm1
datadir <- "../data.chm"
sample_name <- "chm1"
ids <- c(
	"16c58fbcc5633564b10ebe8f78d87883",
	"40c68f29b6d7cb2358f31a7073250406",
	"43a13d07730deb934e9fc01e3b3cd26f",
	"8dcad8fe04f4ebc0ad3254ab4420cdc8",
	"acd889cc16741fb0fba62faa4f7005f3",
	"b1112f1c3cbd28c464f58fc5c5c02f9b",
	"9d134f160ac68c0445002fbb78db4a5e")
truth_id <- "00000000000000000000000000000001"
truth_name <- "Huddleston et al"
generate_figures(datadir, sample_name, ids, truth_id, truth_name)
