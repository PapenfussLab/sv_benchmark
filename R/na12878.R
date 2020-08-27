source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)
library(RColorBrewer)

maxgap <- 180
sizemargin <- 0.50
ignore.strand <- TRUE
minsize <- 51

#####################################
# hg19 blacklist
encodeblacklist <- import(paste0(rootdir, "scripts/input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed"))

#####################################
# na12878 Moleculo/PacBio truth
if (!exists("longsplitreads") || !exists("longspanningreads")) {
	.loadlongreaddelbed <- function(filename) {
		lrbed <- import.bed(con=filename)
		seqlevelsStyle(lrbed) <- "UCSC"
		gr <- GRanges(
			seqnames=c(seqnames(lrbed),seqnames(lrbed)),
			ranges=c(IRanges(start=start(lrbed)-1, width=1), IRanges(start=end(lrbed)+1, width=1)),
			strand=c(rep("+", length(lrbed)), rep("-", length(lrbed))))
		names(gr) <- c(paste0("record", seq_along(lrbed), "_", filename, "/1"), paste0("record", seq_along(lrbed), "_", filename, "/2"))
		gr$partner <- c(paste0("record", seq_along(lrbed), "_", filename, "/2"), paste0("record", seq_along(lrbed), "_", filename, "/1"))
		gr$insLen <- 0
		gr$QUAL <- 0
		gr$FILTER <- "."
		return(gr)
	}
	longsplitreads <- c(
		.loadlongreaddelbed("~/na12878/longread_nov2016/chemistry1.sorted.bam.sr.bam.sr.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/chemistry_2_picard.bam.sr.bam.sr.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/chemistry_3_picard.bam.sr.bam.sr.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sr.bam.sr.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/NA12878.moleculo.bwa-mem.20140110.bam.sr.bam.sr.bed"))
	longspanningreads <- c(
		.loadlongreaddelbed("~/na12878/longread_nov2016/chemistry1.sorted.bam.sp.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/chemistry_2_picard.bam.sp.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/chemistry_3_picard.bam.sp.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sp.bed"),
		.loadlongreaddelbed("~/na12878/longread_nov2016/NA12878.moleculo.bwa-mem.20140110.bam.sp.bed"))
	longsplitreads <- longsplitreads[.distance(longsplitreads, partner(longsplitreads))$max >= minsize * sizemargin]
	longspanningreads <- longspanningreads[.distance(longspanningreads, partner(longspanningreads))$max >= minsize * sizemargin]
}
rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/i/", "~/i/")
#####################################
# Load VCFs
vcfs <- NULL
metadata <- LoadMetadata(paste0(rootdir, "data.na12878"))
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.na12878"), metadata=metadata, existingList=vcfs)
vcfs <- sapply(names(vcfs), function(id) {
	write(paste0("Cleaning ", id), stderr())
	gr <- vcfs[[id]]
	if (length(gr) == 0) return(gr)
	seqlevelsStyle(gr) <- "UCSC"
	# only looking at intrachromosomal events at least 50bp in size
	gr <- gr[!is.na(gr$svLen) & abs(gr$svLen) >= minsize,]
	if (any(is.na(gr$svtype))) {
		warning("NA svtype found - ignoring")
		gr <- gr[!is.na(gr$svtype),]
	}
	# deletion-like events
	gr <- gr[gr$svtype == "DEL" | (gr$svtype == "BND" & strand(gr) != strand(partner(gr)) & strand(gr) == ifelse(start(gr) < start(partner(gr)), "+", "-")),]
	# on primary chromosomes not overlapping blacklist
	gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")),]
	gr <- gr[!overlapsAny(gr, encodeblacklist),]
	gr <- gr[gr$partner %in% names(gr),]
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)

#####################################
# Mills
calls_default <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, truthgr=vcfs[["00000000000000000000000000000002"]])
calls_default$calls$CallSet <- PASS_CALLS
calls_default$truth$CallSet <- PASS_CALLS

calls_all <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, truthgr=vcfs[["00000000000000000000000000000002"]])
calls_all$calls$CallSet <- ALL_CALLS
calls_all$truth$CallSet <- ALL_CALLS
mcalls <- rbind(calls_default$calls, calls_all$calls)

# Duplicates are considered FPs
mcalls$fp[mcalls$duptp] <- TRUE
mcalls$tp[mcalls$duptp] <- FALSE

#####################################
# Long read
lrcalls <- bind_rows(lapply(names(vcfs)[names(vcfs) %in% (metadata %>% filter(!is.na(CX_CALLER)))$Id], function(id) {
	write(paste0("Processing ", id), stderr())
	callgr <- vcfs[[id]]
	result <- data.frame(
		sourceId=callgr$sourceId,
		QUAL=callgr$QUAL,
		FILTER=callgr$FILTER,
		svLen=callgr$svLen) %>% mutate(
		Id=id,
		srhits=0,
		sphits=0)

	.hitCounts <- function(truthgr) {
		hitscounts <- rep(0, length(callgr))
		hits <- findMatchingBreakpoints(callgr, truthgr, maxgap=maxgap, ignore.strand=FALSE, sizemargin=sizemargin)
		hits$QUAL <- callgr$QUAL[hits$queryHits]
		# assign supporting evidence to the call with the highest QUAL
		hits <- hits %>%
			arrange(desc(QUAL), queryHits) %>%
			distinct(subjectHits, .keep_all = TRUE) %>%
			group_by(queryHits) %>%
			summarise(n=n())
		hitscounts[hits$queryHits] <- hits$n
		return(hitscounts)
	}
	result$srhits <- .hitCounts(longsplitreads)
	result$sphits <- .hitCounts(longspanningreads)
	return(result)
}))
lrcalls <- lrcalls %>%
	mutate(tp=srhits>=3 | sphits>=7, fp=!tp, CallSet=ALL_CALLS)
lrcalls <- rbind(lrcalls, lrcalls %>%
	filter(FILTER %in% c(".", "PASS")) %>%
	mutate(CallSet=PASS_CALLS))

#####################################
# Plots
.mostSensitivePerCaller <- function(calls, callers) {
	calls %>%
		dplyr::select(Id, CallSet, tp) %>%
		group_by(Id, CallSet) %>%
		summarise(tp=sum(tp)) %>%
		ungroup() %>%
		arrange(desc(tp)) %>%
		left_join(metadata) %>%
		distinct(CallSet, StripCallerVersion(CX_CALLER), .keep_all = TRUE) %>%
		filter(is.null(callers) | StripCallerVersion(CX_CALLER) %in% callers) %>%
		dplyr::select(Id, CallSet)
}
.plotgraphs <- function(calls, label, callers=c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "manta", "cortex", "hydra")) {
	# use aligner with best sensitivity
	sensAligner <- .mostSensitivePerCaller(calls, callers)
	roc <- calls %>%
		dplyr::select(Id, CallSet, QUAL, tp, fp) %>%
		# force a (0,0) point for all callers
		rbind(calls %>% select(Id, CallSet) %>% distinct(Id, CallSet) %>% mutate(QUAL=2 * max(calls$QUAL), tp=0, fp=0)) %>%
			filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
			arrange(desc(QUAL)) %>%
			group_by(Id, CallSet) %>%
			mutate(tp=cumsum(tp), fp=cumsum(fp)) %>%
			ungroup() %>%
			group_by(Id, CallSet, QUAL) %>%
			summarise(tp=max(tp), fp=max(fp)) %>%
			ungroup() %>%
			mutate(precision=tp / (tp + fp), fdr=1-precision) %>%
			left_join(metadata) %>%
			mutate(caller=StripCallerVersion(CX_CALLER), CallSet=relevel(factor(CallSet), PASS_CALLS))
	write.csv(roc, paste0("na12878_roc_", label, "_error_", maxgap, "bp_", sizemargin, "x", ".csv"))
	ggplot(roc) +
		aes(group=paste(Id, CallSet), y=tp/2, x=fp/2, linetype=CallSet, color=caller) +
		geom_line() +
		coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) +
		scale_colour_brewer(palette="Set3") +
	scale_linetype_manual(values=c("solid", "dotted")) +
		labs(y="True Positives", x="False Positives")
	saveplot(paste0("na12878_tp_fp_", label, "_error_", maxgap, "bp_", sizemargin, "x"), width=7, height=5)
	ggplot(rbind(roc, roc %>%
			# add in horizontal line to the y axis
			filter(tp > 0) %>%
			arrange(tp) %>%
			distinct(caller, CallSet, .keep_all = TRUE) %>%
			mutate(tp=0))) +
		aes(group=paste(Id, CallSet), x=tp/2, y=precision, linetype=CallSet, color=caller) +
		geom_line() +
		coord_cartesian(xlim=c(0, 4000)) +
		scale_color_manual(values=c(brewer.pal("Set1", n=9), "#000000")) +
		#scale_colour_brewer(palette=colorRampPalette("Set1")) +
	scale_linetype_manual(values=c("solid", "dotted")) +
		labs(x="True positives", y="Precision")
	saveplot(paste0("na12878_prec_", label, "_error_", maxgap, "bp_", sizemargin, "x"), width=7, height=5)

	ggplot(rbind(roc, roc %>%
								 # add in horizontal line to the y axis
								 filter(tp > 0) %>%
								 arrange(tp) %>%
								 distinct(caller, CallSet, .keep_all = TRUE) %>%
								 mutate(tp=0)) %>%
								mutate(caller=StripCallerVersion(CX_CALLER)) %>%
								filter(CallSet==PASS_CALLS)) +
		aes(group=paste(Id, CallSet), x=tp/2, y=fdr, color=caller) +
		geom_line() +
		coord_cartesian(xlim=c(0, 3000), ylim=c(0, 0.25)) +
		scale_color_manual(values=c(brewer.pal("Set1", n=9), "#000000")) +
		labs(x="True positives", y="False Discovery Rate (FDR)") +
		theme(plot.margin=unit(c(0,0,0,0), "cm"))
	saveplot(paste0("na12878_fdr_", label, "_error_", maxgap, "bp_", sizemargin, "x"), width=7, height=5)
}

ggplot(lrcalls %>%
			 	filter(tp) %>%
			 	# most sensitive
			 	inner_join(lrcalls %>%
			 						 	group_by(Id, CallSet) %>%
			 						 	summarise(tp=sum(tp)) %>%
			 						 	ungroup() %>%
			 						 	arrange(desc(tp)) %>%
			 						 	left_join(metadata) %>%
			 						 	distinct(CallSet, StripCallerVersion(CX_CALLER), .keep_all = TRUE) %>%
			 						 	dplyr::select(Id, CallSet)
			 						 	) %>%
			 	left_join(metadata) %>%
			 	mutate(caller=StripCallerVersion(CX_CALLER))) +
	aes(group=paste(Id, CallSet), x=abs(svLen), color=caller) +
	scale_x_log10() +
	geom_density() +
	facet_wrap( ~ CallSet)


.plotgraphs(mcalls, "Mills")
.plotgraphs(lrcalls, "pacbiomoleculo")
.plotgraphs(lrcalls, "_all_pacbiomoleculo", callers=fulldatacallers)

dterror <- mcalls %>%
	dplyr::select(Id, CallSet, bperror) %>%
	left_join(metadata) %>%
	mutate(caller=StripCallerVersion(CX_CALLER))
ggplot(dterror) +
	aes(group=paste(Id, CallSet), x=bperror, color=caller, linetype=CallSet) +
	geom_density() +
	scale_linetype_manual(values=c("solid", "dotted")) +
	scale_colour_brewer(palette="Set3") +
	facet_wrap(~caller) +
	labs(title="Error margin", y="density", x="bp error")




