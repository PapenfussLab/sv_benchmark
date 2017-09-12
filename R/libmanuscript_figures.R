source("global.R")
source("sv_benchmark.R")
source("shinyCache2.R")

library(tidyverse)
library(ggplot2)
library(cowplot)
library(colorspace)
library(grid)
library(gridExtra)

library(binom)

## Main figure-generating function ###################################

generate_figures <- function(datadir, sample_name, ids, truth_id, truth_name, grtransformName, longreadbedpedir=NULL, allow_missing_callers=FALSE) {
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

	if (!is.null(longreadbedpedir)) {
		callgr$longreadhits <- 0
		longreadgr <- .CachedLoadTruthBedpe(longreadbedpedir, longread_minMapq)
		lrhits <- findBreakpointOverlaps(callgr, longreadgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin)
		lrhit_summary <- lrhits %>% group_by(queryHits) %>%
			summarise(hitCount=n())
		callgr$longreadhits[lrhit_summary$queryHits] <- lrhit_summary$hitCount
	}

	callgr$truthQUAL <- mcols(callgr)[,paste0("fId",truth_id)]

	callgr$caller_hits_ex_truth <- rowSums(as.matrix(as.data.frame(
		mcols(callgr)[,str_detect(names(mcols(callgr)), "^fId[a-f0-9]+") & !(names(mcols(callgr)) %in% c(paste0("fId", truth_id)))])) != -1)
	callgr$simpleEvent <- simpleEventType(callgr)

	missing_callers <- fulldatacallers[!(fulldatacallers %in% StripCallerVersion((metadata %>% filter(Id %in% callgr$Id))$CX_CALLER))]
	if (!allow_missing_callers && length(missing_callers) > 0) {
		stop(paste("Missing VCF records for ", missing_callers))
	}

	# Figure 1
	plot_overall_roc <- overall_roc_plot(callgr, metadata, truth_id, truth_name)
	saveplot(paste0(fileprefix, "_figure1_roc"), plot=plot_overall_roc, height=6, width=7)

	# Figure 2: simulation
	# done by precache.R

	# Figure 3
	plot3 <- fig_3_grob(callgr, metadata, truth_id, truth_name)
	saveplot(paste0(fileprefix, "_figure3_common_calls"), plot=plot3, height=12, width=14)

	# Figure 5
	plot5 <- fig_5_grob(ids, callgr, metadata)
	saveplot(paste0(fileprefix, "_figure5_qual_bins"), plot=plot5, height=12, width=14)


	# Figure 4:
	plot4 <- roc_by_plots_grob(callgr, metadata, truth_id)
	saveplot(paste0(fileprefix, "_figure4_roc_by"), plot=plot4, height=12, width=14)

	# Supp figures
	plot_dup <- duplicates_ggplot(callgr, truth_id, truth_name, metadata)
	saveplot(paste0(fileprefix, "_Supp_duplicate_call_rate"), plot=plot_dup, height=6, width=7)

	# TODO: call error margin
}

## ROC by ... utilities ##############################################

rocby <- function(callgr, ..., truth_id, rocSlicePoints=100, ignore.duplicates=TRUE) {
	groupingCols <- quos(...)
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
			filter(Id == truth_id & CallSet == ALL_CALLS) %>%
			dplyr::select(!!!groupingCols, dplyr::matches("f?Id.+", ignore.case=FALSE)) %>%
			gather(key="Id_CallSet", value="QUAL", dplyr::matches("f?Id.+", ignore.case=FALSE)) %>%
			mutate(
				fp=0,
				tp=QUAL >= 0,
				fn=QUAL < 0,
				CallSet=ifelse(str_detect(Id_CallSet, "^fId"), ALL_CALLS, PASS_CALLS),
				Id=str_replace(Id_CallSet, "f?Id", "")) %>%
			dplyr::select(-Id_CallSet) %>%
			filter(!fn)
		rocdf <- callgr %>%
			as.data.frame() %>%
			mutate(fp=1, tp=0, fn=0) %>%
			filter(Id != truth_id & truthQUAL < 0) %>%
			dplyr::select(Id, CallSet, !!!groupingCols, QUAL, fp, tp, fn) %>%
			rbind(truthhitsdf)
	}
	return (rocdf %>%
				filter(Id != truth_id) %>%
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

## Plot utilities ####################################################

qual_or_read_count <- function(caller_name) {
	ifelse(
		caller_name %in% c("socrates", "delly", "crest", "pindel", "lumpy", "cortex"),
		"read count",
		"quality score")
	# ... or log quality score! Or log read count ?!?!?!
}

metadata_annotate <- function(df, metadata) {
	df %>%
		left_join(metadata) %>%
		mutate(caller_name = StripCallerVersion(CX_CALLER))
}

n_callers_palette <- function(caller_count, hue = 235) {
	c("black",
	  sequential_hcl(
		  caller_count - 2,
		  h = hue, c. = c(30, 20), l = c(40, 100)),
	  "white") %>%
		rev()
}

caller_colour_scheme <-
	scale_color_manual(
		values = c("#396AB1", "#DA7C30", "#3E9651", "#CC2529",
				   "#535154", "#6B4C9A", "#922428", "#948B3D",
		rep("black", 100)))

#rocby_theme <-
#	cowplot::theme_cowplot() +
#	background_grid("xy", "none", colour.major = "grey70") %+replace%
#	theme(panel.border = element_rect(color = "grey70", linetype = 1, size = 0.2))

roc_common <- function(gg) {
	return(gg +
		aes(y = precision,
			x = tp,
			colour = caller_name,
			linetype = CallSet) +
		geom_line() +
		caller_colour_scheme +
		coord_cartesian(ylim = c(0,1)) +
		scale_y_continuous(labels = scales::percent) +
		theme_cowplot() +
		background_grid("y", "none") +
		labs(
			color = "caller",
			linetype = "call set",
			x = "# true positives",
			y = "precision")
		# Add endpoint showing properties of call set.
	)
}

## Overall ROC ########################################################

overall_roc_plot <- function(callgr, metadata, truth_id, truth_name) {
	plot_out <-
		rocby(callgr, truth_id = truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata) %>%
		ggplot() %>%
		roc_common() +
		labs(title = "Overall Precision-recall")
	return(plot_out)
}

## ROC by stuff ######################################################

roc_by_plots_grob <- function(callgr, metadata, truth_id, genome=str_extract(metadata$CX_REFERENCE, "hg[0-9]+")[1]) {
	callgr$snp50bpbin <- cut(callgr$snp50bp, breaks=c(0, 1, 2, 3, 4, 5, 1000), labels=c(0,1,2,3,4,"5+"), right=FALSE)
	flanking_snvs_rocplot <-
		rocby(callgr, snp50bpbin, truth_id=truth_id) %>%
		metadata_annotate(metadata) %>%
		ggplot() %>%
		roc_common() +
		facet_grid(. ~ snp50bpbin, scales="free") +
		labs(title="Precision-recall by flanking SNV/indels")

	callgr$eventSizeBin <- cut(abs(callgr$svLen), breaks=c(0, 100, 200, 300, 500, 1000, 1000000000), labels=c("50-99", "100-199", "200-299", "300-499", "500-999", "1000+"), right=FALSE)
	eventsize_rocplot <-
		# Is this plot misleading?
		# The number of TP should differ by category -- e.g. "free_x" -- ???.
		rocby(callgr, simpleEvent, eventSizeBin, truth_id=truth_id) %>%
		metadata_annotate(metadata) %>%
		ggplot() %>%
		roc_common() +
		facet_grid(
		    # raw data stratified by simpleEvent
		    . ~ eventSizeBin, scales="free") +
		labs(title="Precision-recall by event size")

	# Unmerged categories; no tandem repeat annotations
	repeatmasker_rocplot <-
		rocby(callgr, repeatClass, simpleEvent, truth_id=truth_id) %>%
		metadata_annotate(metadata) %>%
		ggplot() %>%
		roc_common() +
		facet_wrap(simpleEvent ~ repeatClass, scales="free") +
		labs(title="Precision-recall by RepeatMasker annotation")


	callgr$trf <- overlapsAny(callgr, grtrf[[genome]], type="any")
	callgr$repeatAnn <- ifelse(callgr$repeatClass %in% c("", "DNA", "LINE", "LTR", "SINE", "Other", "Low_complexity", "Simple_repeat"), callgr$repeatClass, "Other")
	callgr$repeatAnn <- ifelse(callgr$repeatAnn == "" & callgr$trf, "TRF", callgr$repeatAnn)
	callgr$repeatAnn <- ifelse(callgr$repeatAnn %in% c("TRF", "Simple_repeat"), "Simple/Tandem", callgr$repeatAnn)
	callgr$repeatAnn <- ifelse(callgr$repeatAnn == "", "No repeat", callgr$repeatAnn)
	callgr$repeatAnn <- str_replace(callgr$repeatAnn, "_", " ")
	callgr$repeatAnn <- relevel(factor(callgr$repeatAnn), ref = "No repeat")
	repeat_rocplot <-
		rocby(callgr, repeatAnn, truth_id=truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata) %>%
		ggplot() %>%
		roc_common() +
		facet_wrap( ~ repeatAnn, scales="free", nrow = 2) +
		labs(title="Precision-recall by presence of repeats at breakpoint")

	# Merge plots

	eventsize_grob <- ggplotGrob(eventsize_rocplot + theme(legend.position = "none"))
	flanking_snvs_grob <- ggplotGrob(flanking_snvs_rocplot + theme(legend.position = "none"))
	repeat_grob <- ggplotGrob(repeat_rocplot)
	fig4_grob <- grid.arrange(
	    grobs = list(eventsize_grob, flanking_snvs_grob, repeat_grob),
	    layout_matrix = rbind(
	        c(1, 1),
	        c(2, 2),
	        c(3, 4)),
	    widths = c(1.50, 0.10),
	    heights = c(1, 1, 1.75))

	return(fig4_grob)

	}


## Figure 3 ##########################################################

shared_tp_calls_plot <- function(callgr, metadata, truth_id, truth_name) {

	truth_hits_df <-
		callgr %>%
		as.data.frame() %>% as.tbl() %>%
		filter(truthQUAL >= 0 | Id == truth_id) %>%
		# Remove duplicates (filtered vs. unfiltered meaningless for truth)
		# Remove non-filtered subject columns
		dplyr::select(Id, CallSet, caller_hits_ex_truth) %>%
		mutate(is_truth = Id==truth_id) %>%
		mutate(truth_factor=factor(1 - is_truth)) %>%
		metadata_annotate(metadata)

	summary_df <- truth_hits_df %>%
		group_by(Id, CallSet, truth_factor) %>%
		summarise(total=n())

	n_callers_plus_truth <-
		length(unique(truth_hits_df$caller_name))

	plot_out <-
		truth_hits_df %>%
		ggplot(aes(x = CallSet)) +
		facet_grid(
			truth_factor + Id ~ .,
			labeller = labeller(
				truth_factor=function(s) {rep("", length(s))},
				Id=function(s) {ifelse(s==truth_id, truth_name,
									   as.character((data.frame(Id=s) %>%
									   metadata_annotate(metadata))$caller_name))}),
			switch = "y") +
		scale_y_continuous(expand = c(0,0)) +
		# scale_x_discrete(labels = )
		geom_bar(aes(fill = interaction(factor(caller_hits_ex_truth), CallSet)),
			size = 0.5, color = "black", width = 1) +
		geom_text(aes(label=total), data=summary_df,
				  hjust=1,
				  y=max(summary_df$total),
				  color="grey50") +
		scale_fill_manual(
			values = c(n_callers_palette(n_callers_plus_truth, 235), n_callers_palette(n_callers_plus_truth, 90)),
			name = "# callers\nsharing",
			labels = rep("", 2 * n_callers_plus_truth)) +
		xlab("") +
		coord_flip() +
		ggtitle("Sharing and distribution of true positive calls") +
		theme_cowplot() +
		# Blanks out "is_truth" panels
		theme(strip.text.y = element_text(angle = 180, hjust = 1),
			  strip.background = element_blank(),
			  axis.text.y = element_blank(),
			  axis.ticks.y = element_blank(),
			  axis.line.y = element_blank())

	return(plot_out)
}

shared_fp_calls_plot <- function(callgr, truth_id, metadata) {

	false_positive_plot_df <-
		callgr[callgr$Id != truth_id] %>%
		as.data.frame() %>% as.tbl() %>%
		filter(truthQUAL < 0) %>%
		# remove duplicate fp calls
		dplyr::select(Id, CallSet, caller_hits_ex_truth, dplyr::matches("f?Id.", ignore.case = FALSE)) %>%
		gather(key = subject_Id_CallSet, value = QUAL, dplyr::matches("f?Id.", ignore.case = FALSE)) %>%
		mutate(subject_Id=str_replace(subject_Id_CallSet, "f?Id", ""),
			   subject_CallSet=ifelse(str_detect(subject_Id_CallSet, "^f"),
									  ALL_CALLS,
									  PASS_CALLS)) %>%
		filter(Id == subject_Id, CallSet == subject_CallSet, QUAL != -2) %>%
		metadata_annotate(metadata)

	n_callers_plus_truth <-
		length(unique(false_positive_plot_df$caller_name))

	ymax_shared <- (false_positive_plot_df %>%
		filter(caller_hits_ex_truth > 1) %>%
		group_by(Id, CallSet) %>%
		summarize(n=n()) %>%
		group_by() %>%
		summarise(n=max(n)))$n

	plot_out <-
		false_positive_plot_df %>%
		ggplot(aes(x = CallSet)) +
		facet_grid(
			Id ~ .,
			labeller = labeller(
				Id=function(s) {
					ifelse(s==truth_id, truth_name,
						as.character((data.frame(Id=s) %>%
							metadata_annotate(metadata))$caller_name))}),
			switch = "y") +
		scale_y_continuous(expand = c(0,0)) +
		# scale_x_discrete(labels = )
		geom_bar(aes(fill = interaction(factor(caller_hits_ex_truth), CallSet)),
				 size = 0.5, color = "black", width = 1) +
		geom_text(aes(label=total), data=false_positive_plot_df %>%
					  group_by(Id, CallSet) %>%
					  summarise(total=n()),
				  y=ymax_shared * 1.1,
				  hjust=1,
				  color="grey50") +
		scale_fill_manual(
			values = c(
				# exclude zero -> white
				# I'm not sure why there needs to be a +1 here
				n_callers_palette(n_callers_plus_truth + 1, 235)[1:n_callers_plus_truth + 1],
				n_callers_palette(n_callers_plus_truth + 1, 90)[1:n_callers_plus_truth + 1]),
			name = "# callers\nsharing",
			labels = rep("", 2 * n_callers_plus_truth)) +
		xlab("") +
		coord_flip(ylim=c(0, 1.1 * ymax_shared)) +
		ggtitle("Sharing and distribution of false positive calls") +
		theme_cowplot() +
		# Blanks out "is_truth" panels
		theme(strip.text.y = element_text(angle = 180, hjust = 1),
			  strip.background = element_blank(),
			  axis.text.y = element_blank(),
			  axis.ticks.y = element_blank(),
			  axis.line.y = element_blank())
	return(plot_out)
}

prec_recall_by_shared_plot <- function(callgr, metadata, truth_id, truth_name) {

	plot_out <-
		rocby(callgr, caller_hits_ex_truth, truth_id = truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata) %>%
		ggplot() %>%
		roc_common() +
		aes(y = precision,
			x = tp,
			colour = caller_name,
			linetype = CallSet) +
		geom_line() +
		facet_wrap(
			~ factor(caller_hits_ex_truth, levels = max(caller_hits_ex_truth):1),
			scale = "free",
			nrow = 2) +
		labs(
			title = "Precision-recall by number of callers sharing call",
			color = "caller",
			linetype = "call set",
			x = "# true positives",
			y = "precision") +
		caller_colour_scheme +
		coord_cartesian(ylim = c(0,1)) +
		scale_y_continuous(labels = scales::percent) +
		theme_cowplot() +
		background_grid("y", "none")

	return(plot_out)
}

fig_3_grob <- function(callgr, metadata, truth_id, truth_name) {

	prec_recall_by_shared_grob <-
		prec_recall_by_shared_plot(
			callgr, metadata, truth_id, truth_name) %>%
		ggplotGrob()

	shared_tp_calls_grob <-
		shared_tp_calls_plot(
			callgr, metadata, truth_id, truth_name) %>%
		ggplotGrob()

	shared_fp_calls_grob <-
		shared_fp_calls_plot(callgr, truth_id, metadata) %>%
		ggplotGrob()

	# shared_calls_grob <-
	#	 rbind(shared_tp_calls_grob, shared_fp_calls_grob)

	layout_matrix <- cbind(c(1,3), c(2,3))

	fig_3_grob <-
		grid.arrange(
			grobs = list(
				shared_tp_calls_grob,
				shared_fp_calls_grob,
				prec_recall_by_shared_grob),
			layout_matrix = layout_matrix)

	return(fig_3_grob)
}

## Figure 5 ##########################################################

get_binned_qual_data <- function(callgr, bin_by, truth_id, bin_count = 100, ci_level = 0.95) {

	df <- data.frame(
		Id = callgr$Id,
		q = callgr$QUAL,
		logq = log10(callgr$QUAL + 1),
		tp = callgr$truthQUAL != -1) %>%
		group_by(Id) %>%
		mutate(
			qbin = as.numeric(as.character(
				cut(q,
					breaks = (0:bin_count)/bin_count * max(q),
					labels = (0:(bin_count - 1))/bin_count * max(q)))),
			logqbin = as.numeric(as.character(
				cut(logq,
					breaks = (0:bin_count)/bin_count * max(logq),
					labels = (0:(bin_count - 1))/bin_count * max(logq))))) %>%
		group_by(Id, .dots = bin_by) %>%
		summarise(
			qmean = mean(q),
			logqmean = mean(logq),
			n = n(),
			prec = sum(tp)/n(),
			prec_lower = binom.confint(sum(tp), n(), ci_level, methods = "exact")$lower,
			prec_upper = binom.confint(sum(tp), n(), ci_level, methods = "exact")$upper)

	return(df)
}


ci_plot <- function(test_id, test_df, qual_column) {

	caller_name <- (metadata %>% filter(Id == test_id))$CX_CALLER

	ci_ggplot <-
		test_df %>%
		ggplot() +
		aes_string(x = qual_column) +
		aes(ymin = prec_lower,
			y = prec,
			ymax = prec_upper) +
		coord_cartesian(ylim = c(0,1)) +
		geom_linerange(
			aes(alpha = factor(n <= 10, levels = c(FALSE, TRUE))),
			color = "grey70") +
		scale_alpha_manual(values = c(1, 0)) +
		geom_point(aes(color = factor(n <= 10, levels = c(FALSE, TRUE)))) +
		scale_color_manual(values = c("#396AB1", "grey70")) +
		cowplot::theme_cowplot() +
		cowplot::background_grid(major = "xy", minor = "none") +
		scale_y_continuous(labels = scales::percent) +
		theme(
			legend.position = "none",
			axis.text.x = element_blank(),
			axis.title.x = element_blank(),
			axis.ticks.x = element_blank(),
			plot.margin = unit(c(1, 1, -0.2, 1), "lines")) +
		ggtitle(StripCallerVersion(caller_name)) +
		ylab("precision")

	return(ci_ggplot)
}

flipped_hist_plot <- function(test_df, qual_column) {

	hist_ggplot <-
		test_df %>%
		ggplot() +
		aes_string(
			x = qual_column,
			xend = qual_column
		) +
		cowplot::theme_cowplot() +
		aes(y = log10(n)) +
		scale_y_reverse(
			expand = c(0,0),
			breaks = 1:6,
			labels = c("10", "", "1k", "", "100k", ""),
			limits = c(7,0)) +
		geom_segment(
			yend = 0,
			color = "grey50") +
		# cowplot::theme_nothing() +
		cowplot::background_grid(major = "xy", minor = "none") +
		theme(
			axis.line = element_blank(),
			axis.ticks = element_blank(),
			plot.margin = unit(c(0, 1, 1, 1), "lines")) +
		ylab("# calls") +
		xlab(str_replace(qual_column, "qmean", " quality score"))

	return(hist_ggplot)
}

stacked_precision_plot <- function(
	test_id, callgr, qual_column = "logqmean") {

	bin_by <- str_replace(qual_column, "mean", "bin")

	test_df <-
		get_binned_qual_data(callgr, bin_by) %>%
		filter(Id == test_id)

	hist_grob <-
		flipped_hist_plot(test_df, qual_column) %>%
		ggplotGrob()

	ci_grob <-
		ci_plot(test_id, test_df, qual_column) %>%
		ggplotGrob()

	combined_grob <-
		gtable_rbind(
			ci_grob, hist_grob,
			size = "first")

	layout_indices <- grep("panel", combined_grob$layout$name)

	panels <- combined_grob$layout$t[layout_indices]

	combined_grob$heights[panels] <-
		combined_grob$heights[panels] * c(2,1)

	return(combined_grob)
}


fig_5_grob <- function(ids, callgr, metadata) {

	all_grobs <-
		map(
			ids,
			(function(caller_id) {
				caller_name <- (metadata %>% filter(Id == caller_id))$CX_CALLER
				stacked_precision_plot(
					caller_id, callgr,
					ifelse(
						StripCallerVersion(caller_name) %in%
							c("manta", "hydra"),
						"qmean",
						"logqmean"))}))

	combined_grob <-
		do.call(
			function(...) grid.arrange(..., ncol = 3),
			all_grobs)

	return(combined_grob)
}

## Duplicates plot ###################################################

duplicates_ggplot <- function(callgr, truth_id, truth_name, metadata) {

	dup_plot_df <- callgr[callgr$Id != truth_id] %>%
		as.data.frame() %>% as.tbl() %>%
		dplyr::select(Id, CallSet, truthQUAL, dplyr::matches("fId")) %>%
		gather(key = subject_Id_CallSet, value = QUAL, dplyr::matches("fId")) %>%
		filter(paste0("fId",Id) == subject_Id_CallSet) %>%
		mutate(
			isTp = truthQUAL != -1,
			isFp = !isTp,
			isDup = QUAL == -2) %>%
		group_by(Id, CallSet) %>%
		summarise(prop_isDup = sum(isDup) / n()) %>%
		ungroup() %>%
		metadata_annotate(metadata)

	dup_plot <- ggplot(dup_plot_df) +
		theme_cowplot() +
		aes(x = caller_name, y = prop_isDup, color = CallSet) +
		geom_errorbar(ymin = 0, ymax = .5, color = "grey50", width = .1) +
		geom_errorbar(ymin = .5, ymax = 1, color = "grey50", width = .1) +
		geom_point(alpha = 0.7, size = 5) +
		scale_color_brewer(palette = "Dark2", name = "") +
		# Actually: percentage of calls overlapping one of higher quality
		labs(x = "", y = "duplicate calls") +
		scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(.05,.05)) +
		theme(axis.line = element_blank(), axis.ticks.y = element_blank(),
			  axis.text.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0))) +
		coord_flip()

	return(dup_plot)
}




















