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

generate_figures_by_eventtype <- function(
	datadir, sample_name, ids, truth_id, truth_name, grtransformName,
	longreadbedpedir=NULL, allow_missing_callers=FALSE) {

	for (eventtype in eventtypes) {
		write(sprintf("Generating figures for event type %s", eventtype), stderr())
		generate_figures(datadir, sample_name, ids, truth_id, truth_name, grtransformName,
										 longreadbedpedir, allow_missing_callers, eventtype)
	}
}

generate_figures <- function(
		datadir, sample_name, ids, truth_id, truth_name, grtransformName,
		longreadbedpedir=NULL, allow_missing_callers=FALSE, eventtype) {

	setCacheRootPath(datadir)
	fileprefix <- str_replace(paste(sample_name, truth_name, ifelse(use_roc_fdr, "fdr", ""), paste0(eventtype, collapse = "_"), sep="_"), "[ /]", "_")
	all_ids <- c(truth_id, ids)
	metadata <- LoadCachedMetadata(datadir)
	metadata <- metadata %>% filter(Id %in% all_ids)
	# force truth
	metadata$CX_REFERENCE_VCF <- list.files(datadir, pattern=paste0("^", truth_id, ".*.vcf$"))


	if (!allow_missing_callers) {
		missing_callers <- fulldatacallers[!(fulldatacallers %in% StripCallerVersion(metadata$CX_CALLER))]
		if (length(missing_callers) > 0) {
			stop(paste("Missing metadata for ", missing_callers))
		}
		vcffiles <- list.files(datadir, pattern="*.vcf")
		missing_vcf <- sapply(ids, function(id) !any(str_detect(vcffiles, id)))
		if (any(missing_vcf)) {
			stop(paste("Missing vcf for ", (data.frame(Id=ids[missing_vcf]) %>% metadata_annotate(metadata))$caller))
		}
	}

	callgr <- .CachedLoadCallMatrixForIds(
		datadir=datadir,
		metadata=metadata,
		ids=all_ids,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand,
		grtransform=lroptions$grtransform[[grtransformName]],
		grtransformName=grtransformName,
		nominalPosition=nominalPosition,
		eventtype=eventtype)
	write(sprintf("callgr generated."), stderr())

	if (length(callgr) == 0) {
		write(paste("No calls for", sample_name, datadir, eventtype), stderr())
		return(NULL)
	}

	callgr$longreadhits <- -1
	if (!is.null(longreadbedpedir)) {
		callgr$longreadhits <- 0
		longreadgr <- .CachedLoadTruthBedpe(longreadbedpedir, longread_minMapq)
		lrhits <- findBreakpointOverlaps(callgr, longreadgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin)
		lrhit_summary <- lrhits %>% group_by(queryHits) %>%
			summarise(hitCount=n())
		callgr$longreadhits[lrhit_summary$queryHits] <- lrhit_summary$hitCount
		write(sprintf("long read annotation complete."), stderr())
	}

	## CALLGR ANNOTATIONS ##

	callgr$truthQUAL <- mcols(callgr)[,paste0("fId",truth_id)]
	callgr$selfQUAL <- self_qual(callgr)

	callgr$caller_hits_ex_truth <- rowSums(as.matrix(as.data.frame(
		mcols(callgr)[,str_detect(names(mcols(callgr)), "^fId[a-f0-9]+") & !(names(mcols(callgr)) %in% c(paste0("fId", truth_id)))])) != -1)
	callgr$simpleEvent <- simpleEventType(callgr)

	# Flanking SNV counts
	callgr$snp50bpbin <-
		cut(callgr$snp50bp,
				breaks = c(0, 1, 2, 3, 4, 5, 1000),
				labels = c(0, 1, 2, 3, 4, "5+"),
				right = FALSE)

	# Binned event sizes
	callgr$eventSizeBin <-
		cut(abs(callgr$svLen),
				breaks = c(0, 100, 200, 300, 500, 1000, 1000000000),
				labels = c("50-99", "100-199", "200-299", "300-499", "500-999", "1000+"),
				right = FALSE)

	# Merged repeat classes
	genome_name <- str_extract(metadata$CX_REFERENCE, "hg[0-9]+")[1]

	callgr$trf <- overlapsAny(callgr, grtrf[[genome_name]], type="any")
	callgr$repeatAnn <- ifelse(callgr$repeatClass %in% c("", "DNA", "LINE", "LTR", "SINE", "Other", "Low_complexity", "Simple_repeat"), callgr$repeatClass, "Other")
	callgr$repeatAnn <- ifelse(callgr$repeatAnn == "" & callgr$trf, "TRF", callgr$repeatAnn)
	callgr$repeatAnn <- ifelse(callgr$repeatAnn %in% c("TRF", "Simple_repeat"), "Simple/Tandem", callgr$repeatAnn)
	callgr$repeatAnn <- ifelse(callgr$repeatAnn == "", "No repeat", callgr$repeatAnn)
	callgr$repeatAnn <- str_replace(callgr$repeatAnn, "_", " ")
	callgr$repeatAnn <- relevel(factor(callgr$repeatAnn), ref = "No repeat")

	# This is a great place from which to debug.
  browser()

	### PLOTTING ###

	write(sprintf("Figure 1"), stderr())
	plot_overall_roc <- overall_roc_plot(callgr, metadata, truth_id, truth_name)
	saveplot(paste0(fileprefix, "_figure1_roc"), plot=plot_overall_roc, height=6, width=7)

	# Figure 2: simulation
	# done by precache.R

	write(sprintf("Figure 5"), stderr())
	plot5 <- fig_5_grob(ids, callgr, metadata)
	saveplot(paste0(fileprefix, "_figure5_qual_bins"), plot=plot5, height=12, width=14)

	write(sprintf("Figure 4"), stderr())
	plot4 <- fig_4_grob(callgr, metadata, truth_id)
	saveplot(paste0(fileprefix, "_figure4_roc_by"), plot=plot4, height=12, width=14)

	write(sprintf("Supp ROC by plots"), stderr())
	# Supp figure:
	roc_by_snp_by_repeats_plot <-
		roc_by_flanking_snvs_by_repeats(callgr, metadata, truth_id, genome)
	saveplot(paste0(fileprefix, "_roc_by_flanking_by_repeat"), plot=roc_by_snp_by_repeats_plot, height=16, width=16)

	# Another supp figure (event size caused by STRs?)
	roc_by_event_size_by_repeats_plot <-
		roc_by_event_size_by_repeats(callgr, metadata, truth_id, genome)
	saveplot(paste0(fileprefix, "_roc_by_event_size_by_repeat"), plot=roc_by_event_size_by_repeats_plot, height=16, width=16)

	write(sprintf("Duplicate call rate"), stderr())
	plot_dup <- duplicates_ggplot(callgr, truth_id, truth_name, metadata)
	saveplot(paste0(fileprefix, "_Supp_duplicate_call_rate"), plot=plot_dup, height=4.5, width=7)

	write(sprintf("Figure 3"), stderr())
	plot3 <- fig_3_grob(callgr, metadata, truth_id, truth_name)
	saveplot(paste0(fileprefix, "_figure3_common_calls"), plot=plot3, height=12, width=14)

	write(sprintf("Ensemble"), stderr())
	plot_ensemble <- ensemble_plot_list(callgr, metadata, truth_id, ids)
	# saveplot(paste0(fileprefix, "_ensemble"), plot=plot_ensemble, height=12, width=14)

	# TODO: call error margin

	# TODO: this should be empty
	table((callgr$selfQUAL < callgr$QUAL) & callgr$QUAL >= 0)
}

## ROC by ... utilities ##############################################

rocby <- function(callgr, ..., truth_id, rocSlicePoints=100, ignore.duplicates=TRUE, minlongreadhits=1000000000) {
	groupingCols <- quos(...)
	if (!ignore.duplicates) {
		rocdf <- callgr %>%
			as.data.frame() %>%
			dplyr::select(Id, CallSet, !!!groupingCols, QUAL, truthQUAL, longreadhits) %>%
			mutate(tp=truthQUAL > 0 | longreadhits >= minlongreadhits,
						 fp=!tp) %>%
			select(-longreadhits)
	} else {
		isdup <- callgr$Id != truth_id & # never dedup truth calls
			(callgr$truthQUAL == -2 | # dedup tp hits
			callgr$truthQUAL == -1 & callgr$selfQUAL == -2) # dedup fp misses

		# use the truth attributes for hits
		truthhitsdf <- callgr %>%
			as.data.frame() %>%
			filter(Id == truth_id & CallSet == ALL_CALLS) %>%
			dplyr::select(!!!groupingCols, dplyr::matches("f?Id.+", ignore.case=FALSE)) %>%
			gather(key="Id_CallSet", value="QUAL", dplyr::matches("f?Id.+", ignore.case=FALSE)) %>%
			mutate(
				fp=FALSE,
				tp=QUAL >= 0,
				fn=QUAL < 0,
				CallSet=colname_to_CallSet(Id_CallSet),
				Id=colname_to_Id(Id_CallSet)) %>%
			dplyr::select(-Id_CallSet)

		islrtp_only <- callgr$truthQUAL == -1 & callgr$longreadhits >= minlongreadhits
		longread_truthhitsdf <- callgr[islrtp_only] %>%
			as.data.frame() %>%
			filter(!(Id == truth_id)) %>%
			mutate(fp=FALSE, tp=TRUE, fn=FALSE) %>%
			dplyr::select(Id, CallSet, !!!groupingCols, QUAL, fp, tp, fn)

		rocdf <- callgr[!islrtp_only] %>%
			as.data.frame() %>%
			mutate(fp=TRUE, tp=FALSE, fn=FALSE) %>%
			filter(Id != truth_id, truthQUAL < 0) %>%
			dplyr::select(Id, CallSet, !!!groupingCols, QUAL, fp, tp, fn) %>%
			bind_rows(truthhitsdf, longread_truthhitsdf)
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
				mutate(is_endpoint = tp == max(tp)) %>%
				ungroup() %>%
				#TODO: fn and sens need eventCount using group_by(Id, CallSet, !!!groupingCols)
				mutate(
					precision=tp / (tp + fp),
					fdr=1-precision))
}

self_qual <- function(callgr) {
	qualcols <- names(mcols(callgr)) %>%
		(function(x) {x[str_detect(x, "^.?Id.+")]})
	quals <- as.matrix(mcols(callgr)[,qualcols])
	qual_col_name_matrix <- matrix(rep(qualcols, length(callgr)), ncol=length(qualcols), byrow=TRUE)
	is_self_col <- ifelse(qual_col_name_matrix == IdCallSet_to_colname(callgr$Id, callgr$CallSet), 1, 0)
	return(rowSums(quals * is_self_col))
}

## Plot utilities ####################################################

display_log_qual <- function(caller) {
	!(StripCallerVersion(caller) %in% c("manta", "hydra"))
}

qual_or_read_count <- function(caller) {
	paste0(ifelse(StripCallerVersion(caller) %in% c(
		"socrates", "delly", "crest", "pindel", "lumpy", "cortex"),
		"read count", "quality score"),
				 ifelse(display_log_qual(caller), " + 1", "")
	)
}

n_callers_palette <- function(caller_count, hue = 235, crange = c(30, 20), lrange = c(40, 100)) {
	c("black",
	  sequential_hcl(
		  caller_count - 2,
		  h = hue, c. = crange, l = lrange),
	  "white") %>%
		rev()
}

fixed_caller_metadata <- data_frame(
	caller_name=
		c("pindel",  "manta",   "breakdancer", "hydra",   "gridss",  "socrates", "crest",   "cortex",  "delly",   "lumpy"),
	caller_initial =
		c("P",       "M",       "B",           "H",       "G",       "S",        "C",       "X",       "D",       "L"),
	caller_colour =
		c("#CC2529", "#396AB1", "#3E9651",     "#45b0cd", "#535154", "#6B4C9A",  "#922428", "#948B3D", "#e2a198", "#DA7C30")
)

caller_colour_scheme <-
	scale_colour_manual(
		values = purrr::set_names(
			fixed_caller_metadata$caller_colour,
			fixed_caller_metadata$caller_name),
		na.value = "black"
	)

metadata_annotate <- function(df, metadata) {
	df %>%
		left_join(metadata, by="Id") %>%
		mutate(caller_name = StripCallerVersion(CX_CALLER)) %>%
		left_join(fixed_caller_metadata) %>%
		replace_na(list(
			caller_colour = "black",
			caller_initial = "﹡",
			caller_name = "(no name)"
		)) %>%
		mutate(
			caller_name = factor(caller_name)
		)
}

use_roc_fdr <- FALSE
roc_title <- function() { ifelse(use_roc_fdr, "FDR-recall", "Precision-recall")}

roc_common <- function(df) {
	gg <- ggplot(df) +
		aes(
			y = ifelse(use_roc_fdr, 1 - precision, precision),
			x = tp,
			colour = caller_name,
			linetype = CallSet) +
		geom_line(size = 0.3) +
		# Baubles -- optional?
		geom_point(data = df %>% filter(is_endpoint), size = 6
							 # , colour = "black", fill = "white", stroke = 0.3, shape = 21
							 ) +
		geom_text(
			aes(label = caller_initial),
			data = df %>% filter(is_endpoint),
			# nudge_x = .2, nudge_y = .2,
			fontface = "bold",
			hjust = "center", vjust = "center",
			nudge_y = .001,
			color = "white") +
		caller_colour_scheme +
		coord_cartesian(ylim = c(0,1)) +
		scale_y_continuous(labels = scales::percent) +
		theme_cowplot() +
		background_grid("y", "none") +
		labs(
			color = "caller",
			linetype = "call set",
			x = "# true positives"
			)
		# Add endpoint showing properties of call set.
	if (use_roc_fdr) {
		gg <- gg +
			aes(y = fdr) +
			labs(y = "false discovery rate")
	} else {
		gg <- gg +
			aes(y = precision) +
			labs(y = "precision")
	}
	return(gg)
}

## Overall ROC ########################################################

overall_roc_plot <- function(callgr, metadata, truth_id, truth_name) {
	plot_out <-
		rocby(callgr, truth_id = truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata) %>%
		roc_common() +
		labs(title = roc_title())
	return(plot_out)
}

## ROC by stuff ######################################################

roc_by_flanking_snvs <- function(callgr, metadata, truth_id) {

	flanking_snvs_rocplot <-
		rocby(callgr, snp50bpbin, truth_id = truth_id) %>%
		metadata_annotate(metadata) %>%
		roc_common() +
		facet_grid(. ~ snp50bpbin, scales="free") +
		labs(title=paste(roc_title(), "by flanking SNV/indels"))

	return(flanking_snvs_rocplot)
}

roc_by_eventsize <- function(callgr, metadata, truth_id) {

	eventsize_rocplot <-
		# Is this plot misleading?
		# The number of TP should differ by category -- e.g. "free_x" -- ???.
		rocby(callgr, simpleEvent, eventSizeBin, truth_id = truth_id) %>%
		metadata_annotate(metadata) %>%
		roc_common() +
		facet_grid(
			# raw data stratified by simpleEvent
			. ~ eventSizeBin, scales="free") +
		labs(title=paste(roc_title(), "by event size"))

	return(eventsize_rocplot)
}

# Doesn't appear in figure
roc_by_repeatmasker <- function(callgr, metadata, truth_id) {

	# Unmerged categories; no tandem repeat annotations
	repeatmasker_rocplot <-
		rocby(callgr, repeatClass, simpleEvent, truth_id=truth_id) %>%
		metadata_annotate(metadata) %>%
		roc_common() +
		facet_wrap(simpleEvent ~ repeatClass, scales="free") +
		labs(title=paste(roc_title(), "by RepeatMasker annotation"))

	return(repeatmasker_rocplot)
}

roc_by_repeat_class_merged <- function(callgr, metadata, truth_id, genome) {

	repeat_rocplot <-
		rocby(callgr, repeatAnn, truth_id=truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata) %>%
		roc_common() +
		facet_wrap( ~ repeatAnn, scales="free", nrow = 2) +
		labs(title=paste(roc_title(), "by presence of repeats at breakpoint"))

	return(repeat_rocplot)
}

# ROC by (flanking x repeat)
# Doesn't appear in final plot
roc_by_flanking_snvs_by_repeats <- function(callgr, metadata, truth_id, genome) {

	grouped_plot_df <-
		rocby(callgr, snp50bpbin, repeatAnn, truth_id = truth_id) %>%
		# Need to filter truth_id here, as in roc_by_repeat_class_merged?
		metadata_annotate(metadata) %>%
		group_by(snp50bpbin, repeatAnn) %>% # Maybe this should only be grouped by one of these??
		mutate(
			scaled_tp=tp/max(tp),
			max_tp = max(tp))

	roc_by_flanking_snvs_by_repeats_plot <-
		grouped_plot_df %>%
		ungroup() %>%
		roc_common() +
			aes(x=scaled_tp) +
		facet_grid(repeatAnn ~ snp50bpbin, scales="free") +
		labs(title = paste(roc_title(), "by presence of repeats at breakpoint\nand flanking SNV/indels"),
				 x = "relative sensitivity") +
		geom_text(
			data = grouped_plot_df %>% distinct(snp50bpbin, repeatAnn, max_tp, .keep_all = TRUE),
			aes(label = max_tp),
			x = .5, y = .5, size = 12, color = "grey40", alpha = .6)



	return(roc_by_flanking_snvs_by_repeats_plot)
}

# ROC by (SV len (eventsizebin) x repeat)
# Doesn't appear in final plot
roc_by_event_size_by_repeats <- function(callgr, metadata, truth_id, genome) {

	grouped_plot_df <-
		rocby(callgr, eventSizeBin, repeatAnn, truth_id = truth_id) %>%
		# Need to filter truth_id here, as in roc_by_repeat_class_merged?
		metadata_annotate(metadata) %>%
		group_by(eventSizeBin, repeatAnn) %>% # Maybe this should only be grouped by one of these??
		mutate(
			scaled_tp=tp/max(tp),
			max_tp = max(tp))

	roc_by_event_size_by_repeats_plot <-
		grouped_plot_df %>%
		ungroup() %>%
		roc_common() +
		aes(x=scaled_tp) +
		facet_grid(repeatAnn ~ eventSizeBin, scales="free") +
		labs(title = paste(roc_title(), "by presence of repeats at breakpoint\nand event size"),
				 x = "relative sensitivity") +
		geom_text(
			data = grouped_plot_df %>% distinct(eventSizeBin, repeatAnn, max_tp, .keep_all = TRUE),
			aes(label = max_tp),
			x = .5, y = .5, size = 12, color = "grey40", alpha = .6)

	return(roc_by_event_size_by_repeats_plot)
}

fig_4_grob <- function(callgr, metadata, truth_id,
											 genome = str_extract(metadata$CX_REFERENCE, "hg[0-9]+")[1]) {
	# Merge plots

	eventsize_grob <-
		ggplotGrob(roc_by_eventsize(callgr, metadata, truth_id) +
			theme(legend.position = "none"))

	flanking_snvs_grob <-
		ggplotGrob(roc_by_flanking_snvs(callgr, metadata, truth_id) +
		theme(legend.position = "none")
		)

	repeat_grob <-
		roc_by_repeat_class_merged(callgr, metadata, truth_id, genome) %>%
		ggplotGrob()

	fig4_grob <- grid.arrange(
		grobs = list(eventsize_grob, flanking_snvs_grob, repeat_grob),
		layout_matrix = rbind(
			c(2, 2),
			c(3, 4),
			c(1, 1)),
		widths = c(1.50, 0.10),
		heights = c(1, 1.75, 1))

	return(fig4_grob)
}

## Figure 3 ##########################################################

make_shared_tp_calls_grob <- function(callgr, metadata, truth_id, truth_name) {

	truth_hits_df <-
		callgr %>%
		as.data.frame() %>% as.tbl() %>%
		filter(truthQUAL >= 0 | Id == truth_id) %>%
		# Remove duplicates (filtered vs. unfiltered meaningless for truth)
		# Remove non-filtered subject columns
		dplyr::select(Id, CallSet, caller_hits_ex_truth) %>%
		mutate(is_truth = Id==truth_id) %>%
		metadata_annotate(metadata)

	p_empty <- grid.text("truth_hits_df empty!")

	if (nrow(truth_hits_df) == 0) {
		return(p_empty)
	}

	summary_df <-
		truth_hits_df %>%
		group_by(Id, caller_name, CallSet) %>%
		summarise(total=n()) %>%
		ungroup()

	# Is the "PASS" field distinct for each caller?
	distinct_pass_df <-
		summary_df %>%
		select(Id, caller_name, CallSet, total) %>%
		spread(key = CallSet, value = total) %>%
		mutate(distinct_pass = `All calls` != `PASS only`) %>%
		select(Id, caller_name, distinct_pass)

	summary_df_condensed <-
		summary_df %>%
		left_join(distinct_pass_df) %>%
		filter((CallSet == "All calls" | distinct_pass)) %>%
		mutate(
			caller_name =
				caller_name %>%
				relevel("NA") %>%
				recode(`NA` = paste("(", truth_name, ")", sep = "")))

	n_callers_plus_truth <-
		length(unique(truth_hits_df$caller_name))

	plot_df <-
		truth_hits_df %>%
		left_join(distinct_pass_df) %>%
		mutate(
			caller_name =
				caller_name %>%
				relevel("NA") %>%
				recode(`NA` = paste("(", truth_name, ")", sep = ""))) %>%
		filter(
			# Show only "All calls" unles there's a meaningful "PASS" field
			(CallSet == "All calls" | distinct_pass))

	plot_out <-
		plot_df %>%
		ggplot(aes(x = CallSet)) +
		facet_grid(
			caller_name ~ .,
			switch = "y", scales = "free", space = "free") +
		scale_y_continuous(expand = c(0,0)) +
		# scale_x_discrete(labels = )
		geom_bar(aes(fill = interaction(factor(caller_hits_ex_truth), CallSet)),
			size = 0.5, color = "black", width = 1) +
		geom_text(aes(label = total), data = summary_df_condensed,
				  hjust = 1,
				  y = max(summary_df$total),
				  color = "grey50") +
		scale_fill_manual(
			values = c(
				n_callers_palette(n_callers_plus_truth, 235),
				# we need to avoid "white" because the truth set isn't "PASS"
				n_callers_palette(n_callers_plus_truth, 90)[2:(n_callers_plus_truth)]),
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

	grob_out <- ggplotGrob(plot_out)

	return(grob_out)
}

make_shared_fp_calls_grob <- function(callgr, truth_id, metadata) {

	false_positive_plot_df <-
		callgr[callgr$Id != truth_id] %>%
		as.data.frame() %>% as.tbl() %>%
		filter(truthQUAL < 0, selfQUAL != -2) %>% # remove tp & duplicate fp calls
		dplyr::select(Id, CallSet, caller_hits_ex_truth) %>%
		metadata_annotate(metadata)


	if (nrow(false_positive_plot_df) == 0) {
		return(grid.text("false_positive_plot_df empty!"))
	}

	summary_df <- false_positive_plot_df %>%
		group_by(Id, caller_name, CallSet) %>%
		summarise(total=n())

	# Is the "PASS" field distinct for each caller?
	distinct_pass_df <-
		summary_df %>%
		select(Id, caller_name, CallSet, total) %>%
		spread(key = CallSet, value = total) %>%
		mutate(distinct_pass = `All calls` != `PASS only`) %>%
		select(Id, caller_name, distinct_pass)

	n_callers_plus_truth <-
		length(unique(false_positive_plot_df$caller_name))

	ymax_shared <- (false_positive_plot_df %>%
		filter(caller_hits_ex_truth > 1) %>%
		group_by(Id, CallSet) %>%
		summarize(n=n()) %>%
		group_by() %>%
		summarise(n=max(n)))$n

	plot_df <-
		false_positive_plot_df %>%
		left_join(distinct_pass_df) %>%
		filter(
			# Show only "All calls" unles there's a meaningful "PASS" field
			(CallSet == "All calls" | distinct_pass))

	plot_out <-
		plot_df %>%
		ggplot(aes(x = CallSet)) +
		facet_grid(
			caller_name ~ .,
			switch = "y", scales = "free", space = "free") +
		scale_y_continuous(expand = c(0,0)) +
		# scale_x_discrete(labels = )
		geom_bar(aes(fill = interaction(factor(caller_hits_ex_truth), CallSet)),
				 size = 0.5, color = "black", width = 1) +
		geom_text(
			aes(label = total),
			data = (summary_df %>%
							 left_join(distinct_pass_df) %>%
							 filter((CallSet == "All calls" | distinct_pass))),
		  y = ymax_shared * 1.1,
		  hjust = 1,
		  color = "grey50") +
		scale_fill_manual(
			values = c(
				# exclude zero -> white
				# I'm not sure why there needs to be a +1 here
				n_callers_palette(n_callers_plus_truth + 1, 235)[1:n_callers_plus_truth + 1],
				n_callers_palette(n_callers_plus_truth + 1, 90 )[1:n_callers_plus_truth + 1]),
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

	grob_out <- ggplotGrob(plot_out)

	return(grob_out)
}

prec_recall_by_shared_plot <- function(callgr, metadata, truth_id, truth_name) {

	plot_out <-
		rocby(callgr, caller_hits_ex_truth, truth_id = truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata) %>%
		roc_common() +
		facet_wrap(
			~ factor(caller_hits_ex_truth, levels = max(caller_hits_ex_truth):1),
			scale = "free",
			nrow = 2) +
		labs(
			title = paste(roc_title(), "by number of callers sharing call"))

	return(plot_out)
}

fig_3_grob <- function(callgr, metadata, truth_id, truth_name) {

	prec_recall_by_shared_grob <-
		prec_recall_by_shared_plot(
			callgr, metadata, truth_id, truth_name) %>%
		ggplotGrob()

	shared_tp_calls_grob <-
		make_shared_tp_calls_grob(
			callgr, metadata, truth_id, truth_name)

	shared_fp_calls_grob <-
		make_shared_fp_calls_grob(callgr, truth_id, metadata)

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

get_binned_qual_data <- function(callgr, bin_by, bin_count = 100, ci_level = 0.95) {

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


ci_plot <- function(test_id, test_df, qual_column, metadata) {

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

flipped_hist_plot <- function(test_df, qual_column, caller_name) {

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
		xlab(qual_or_read_count(caller_name))

	if (str_detect(qual_column, "^log")) {
		max_log_qual_level <- max(test_df[[qual_column]])

		all_labels <- sort(c(10**(0:(max_log_qual_level + 1)), # 1, 10, ...
												 10**(0:(max_log_qual_level + 1)) * 3)) # 3, 30, ...

		labels <- all_labels[all_labels < 10**max_log_qual_level]

		if (length(labels) > 7) {
			labels <- (10**(0:(max_log_qual_level + 1)))
			if (length(labels) > 7) {
				labels <- labels[seq(1, length(labels), 2)]
			}
			breaks <- log10(labels)
			labels <- scales::trans_format('log10', scales::math_format())(labels)
		} else {
			breaks <- log10(labels)
		}

		hist_ggplot <- hist_ggplot + scale_x_continuous(breaks = breaks, labels = labels)
	}

	return(hist_ggplot)
}

stacked_precision_plot <- function(
	test_id, callgr, metadata, qual_column = "logqmean", caller_name) {

	bin_by <- str_replace(qual_column, "mean", "bin")

	test_df <-
		get_binned_qual_data(callgr, bin_by) %>%
		filter(Id == test_id)

	if (nrow(test_df) == 0) {
		return(grid.text("test_df empty!"))
	}

	hist_grob <-
		flipped_hist_plot(test_df, qual_column, caller_name) %>%
		ggplotGrob()

	ci_grob <-
		ci_plot(test_id, test_df, qual_column, metadata) %>%
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
					caller_id, callgr, metadata,
					ifelse(display_log_qual(caller_name), "logqmean", "qmean"),
					caller_name
					)}))

	combined_grob <-
		do.call(
			function(...) grid.arrange(..., ncol = 3),
			all_grobs)

	return(combined_grob)
}

## Duplicates plot ###################################################

duplicates_ggplot <- function(callgr, truth_id, truth_name, metadata) {

	# browse()

	dup_plot_df <- callgr[callgr$Id != truth_id] %>%
		as.data.frame() %>% as.tbl() %>%
		dplyr::select(Id, CallSet, truthQUAL, selfQUAL) %>%
		mutate(
			isTp = truthQUAL != -1,
			isFp = !isTp,
			isDup = selfQUAL == -2) %>%
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

## Ensemble calling plot ###################################################

#' @params p number of callers to calculate ensemble for (nCp)
ensemble_plot_list <- function(callgr, metadata, truth_id, ids, p=length(ids), minlongreadhits=1000000000) {

	ensemble_df <- calc_ensemble_performance(callgr, metadata, ids, p, minlongreadhits)
	eventcount <- sum(callgr$Id == truth_id & callgr$CallSet == ALL_CALLS)
	ensemble_df <-
		ensemble_df %>%
		mutate(
			sens = tp / eventcount,
			f1score = 2 * precision * sens / (precision + sens),
			ensemble = str_c(minhits, " of ", ncallers))

	faceted_plot <-
		ggplot(ensemble_df) +
		aes(x=tp, y=precision, colour=CallSet) + #, shape=as.factor(minhits))
		geom_point() +
		coord_cartesian(ylim = c(0,1)) +
		scale_y_continuous(labels = scales::percent) +
		theme_cowplot() +
		background_grid("xy", "none") +
		facet_grid(ncallers ~ minhits) +
		labs(title="Ensemble caller performance for calls requiring x of y callers") +
		theme_cowplot() +
		scale_color_brewer(palette = "Set2") +
		background_grid()

	pareto_frontier_by_ensemble_df <-
		ensemble_df %>%
		group_by(minhits, ncallers) %>%
		arrange(desc(precision)) %>%
		filter(sens == cummax(sens))

	pareto_frontier_df <-
		ensemble_df %>%
		arrange(desc(precision)) %>%
		filter(sens == cummax(sens))

	overall_roc_df <-
		rocby(callgr, truth_id = truth_id) %>%
		filter(Id != truth_id) %>%
		metadata_annotate(metadata)

	bauble_points <-
		geom_point(data = overall_roc_df %>% filter(is_endpoint), size = 6)

	bauble_text <-
		geom_text(
			aes(label = caller_initial),
			data = overall_roc_df %>% filter(is_endpoint),
			# nudge_x = .2, nudge_y = .2,
			fontface = "bold",
			hjust = "center", vjust = "center",
			nudge_y = .001,
			color = "white")

	overall_roc_plot <-
		ggplot(overall_roc_df) +
		aes(
			y = precision,
			x = tp) +
		coord_cartesian(ylim = c(0,1)) +
		scale_y_continuous(labels = scales::percent) +
		theme_cowplot() +
		background_grid("y", "none")

	pareto_frontier_plot_overall <-
		overall_roc_plot +
		geom_point(
			data = ensemble_df %>% filter(ensemble %in% c("1 of 5", "2 of 3", "4 of 5")),
			aes(color = ensemble),
			size = 0.6) +
		geom_line(
			data = pareto_frontier_df) +
		scale_color_brewer(palette = "Set2") +
		bauble_points + bauble_text

	all_ensemble_plot <-
		overall_roc_plot +
		geom_point(
			data = ensemble_df %>% mutate(minhits = pmin(minhits, 3), ncallers = pmin(ncallers, 4), ensemble = str_c(minhits, ifelse(minhits >= 3, "+", ""), " of ", ncallers, ifelse(ncallers >= 4, "+", ""))) %>% filter(ensemble != "1 of 1") %>% mutate(ensemble=ifelse(minhits==1, "1 of n", ensemble)) %>% arrange(desc(minhits)),
			aes(color = ensemble),
			size = 0.6) +
		scale_color_brewer(palette = "Set2") +
		bauble_points + bauble_text

	# ggplot() +
	# 	aes(x=tp, y=precision) +
	# 	geom_point(data=ensemble_df,
	# 		aes(colour=factor(ncallers), shape=factor(minhits))) +
	# 	geom_line(data=rocby(callgr, truth_id=truth_id, minlongreadhits=minlongreadhits) %>%
	# 			metadata_annotate(metadata),
	# 		aes(group=interaction(Id, CallSet)),
	# 		colour="grey50", size=2) +
	# 	scale_shape_manual(values=1:max(ensemble_df$minhits)) +
	# 	labs(title="Ensemble caller performance for calls requiring x of y callers") +
	# 	geom_point(data=ensemble_df %>% filter(
	# 		ncallers == 3,
	# 		minhits %in% c(2,3),
	# 		str_detect(callers, "delly"),
	# 		str_detect(callers, "manta"),
	# 		str_detect(callers, "lumpy"),
	# 		str_detect(callers, "delly")),
	# 		size=3,
	# 		colour="black") +
	# 	geom_text(data=ensemble_df %>% filter(
	# 		ncallers == 3,
	# 		minhits %in% c(2,3),
	# 		str_detect(callers, "delly"),
	# 		str_detect(callers, "manta"),
	# 		str_detect(callers, "lumpy"),
	# 		str_detect(callers, "delly")),
	# 		size=3,
	# 		colour="black",
	# 		aes(label=callers))
#
# 	ggplot() +
# 		aes(x=tp, y=precision) +
# 		geom_point(data=ensemble_df %>% filter(!str_detect(callers, "gridss") & !str_detect(callers, "manta")),
# 							 aes(colour=factor(ncallers), shape=factor(minhits))) +
# 		geom_line(data=rocby(callgr, truth_id=truth_id, minlongreadhits=minlongreadhits) %>%
# 								metadata_annotate(metadata),
# 							aes(group=interaction(Id, CallSet)),
# 							colour="grey50", size=0.5) +
# 		scale_shape_manual(values=1:max(ensemble_df$minhits)) +
# 		labs(title="Ensemble caller performance for calls requiring x of y callers\nEnsembles using GRIDSS or manta excluded.")
#
# 	ggplot() +
# 		aes(x=tp, y=precision) +
# 		#geom_point(data=ensemble_df %>% filter(str_detect(callers, "gridss") | str_detect(callers, "manta")),
# 		#					 aes(shape=factor(minhits)), colour="black") +
# 		geom_point(data=ensemble_df %>% filter(!str_detect(callers, "gridss") & !str_detect(callers, "manta")),
# 							 aes(colour=factor(ncallers), shape=factor(minhits))) +
# 		geom_line(data=rocby(callgr, truth_id=truth_id, minlongreadhits=minlongreadhits) %>%
# 								metadata_annotate(metadata),
# 							aes(group=interaction(Id, CallSet)),
# 							colour="grey50", size=0.5) +
# 		scale_shape_manual(values=1:max(ensemble_df$minhits)) +
# 		labs(title="Ensemble caller performance for calls requiring x of y callers\nEnsembles using GRIDSS or manta excluded.") +
# 		geom_point(data=ensemble_df %>% filter(
# 			ncallers == 3,
# 			minhits %in% c(2,3),
# 			str_detect(callers, "delly"),
# 			str_detect(callers, "manta"),
# 			str_detect(callers, "lumpy"),
# 			str_detect(callers, "delly")),
# 			size=3,
# 			colour="black") +
# 		geom_text(data=ensemble_df %>% filter(
# 			ncallers == 3,
# 			minhits %in% c(2,3),
# 			str_detect(callers, "delly"),
# 			str_detect(callers, "manta"),
# 			str_detect(callers, "lumpy"),
# 			str_detect(callers, "delly")),
# 			size=3,
# 			colour="black",
# 			aes(label=callers))

	return(list(
		pareto_frontier_plot_overall,
		all_ensemble_plot,
		faceted_plot
	))
}

calc_ensemble_performance <- function(callgr, metadata, ids, p, minlongreadhits) {
	allgr <- callgr[callgr$CallSet==ALL_CALLS]
	passgr <- callgr[callgr$CallSet==PASS_CALLS]
	bind_rows(lapply(1:p, function(choosep) {
		id_subset <- combn(ids, choosep)
		write(sprintf("Calculating ensemble calls nC%d", choosep), stderr())
		bind_rows(lapply(1:ncol(id_subset), function(i) {
			ensemble_ids <- id_subset[,i]
			bind_rows(
				.ensemble_performance(allgr, ensemble_ids, ALL_CALLS, minlongreadhits),
				.ensemble_performance(passgr, ensemble_ids, PASS_CALLS, minlongreadhits)
			) %>% mutate(
				callers=paste(StripCallerVersion((data.frame(Id=ensemble_ids) %>% metadata_annotate(metadata))$CX_CALLER), collapse=","))
		}))
	})) %>%
		mutate(
			precision=tp / (tp + fp),
			fdr=1-precision)
}
.ensemble_performance <- function(callgr, ensemble_ids, ensemble_callset, minlongreadhits) {
	# for each call, count the number ensemble callers that called it
	hits <- rowSums(as.matrix(mcols(callgr)[,IdCallSet_to_colname(ensemble_ids, ensemble_callset)]) != -1)
	bind_rows(lapply(1:length(ensemble_ids), function(requiredHits) {
		relevant_call <-
			hits >= requiredHits &
			callgr$Id %in% ensemble_ids &
			callgr$CallSet == ensemble_callset &
			(callgr$truthQUAL >= 0 | callgr$selfQUAL != -2) # exclude self duplicates that aren't TPs
		tp <- callgr$truthQUAL != -1 | callgr$longreadhits >= minlongreadhits
		return(data.frame(
			ncallers=length(ensemble_ids),
			minhits=requiredHits,
			# To count calls we count all the calls across all callers that are in the intersection
			# To prevent overcounting (a call made by 3 callers will be in the call set three times)
			# we divide out by the number of relevant callers making the call (3 callers making the same fp * 1/3 = 1 fp)
			tp=sum(1/hits[relevant_call & tp]),
			fp=sum(1/hits[relevant_call & !tp]),
			CallSet=ensemble_callset))
	}))
}



















