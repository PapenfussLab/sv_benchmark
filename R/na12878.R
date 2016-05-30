source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)

maxgap <- 180
sizemargin <- 0.5
ignore.strand <- TRUE
minsize <- 51

#####################################
# hg19 blacklist
encodeblacklist <- import(paste0(rootdir, "input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed"))

#####################################
# na12878 Moleculo/PacBio truth
if (!exists("longsplitreads") || !exists("longspanningreads")) {
	.loadlongreaddelbed <- function(filename) {
		lrbed <- import.bed(con=paste0(rootdir, filename))
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
	  .loadlongreaddelbed("input.na12878/chemistry1.sorted.bam.sr.bam.sr.bed"),
	  .loadlongreaddelbed("input.na12878/chemistry_2_picard.bam.sr.bam.sr.bed"),
	  .loadlongreaddelbed("input.na12878/chemistry_3_picard.bam.sr.bam.sr.bed"),
	  .loadlongreaddelbed("input.na12878/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sr.bam.sr.bed"),
	  .loadlongreaddelbed("input.na12878/NA12878.moleculo.bwa-mem.20140110.bam.sr.bam.sr.bed"))
	longspanningreads <- c(
	  .loadlongreaddelbed("input.na12878/chemistry1.sorted.bam.sp.bed"),
	  .loadlongreaddelbed("input.na12878/chemistry_2_picard.bam.sp.bed"),
	  .loadlongreaddelbed("input.na12878/chemistry_3_picard.bam.sp.bed"),
	  .loadlongreaddelbed("input.na12878/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sp.bed"),
	  .loadlongreaddelbed("input.na12878/NA12878.moleculo.bwa-mem.20140110.bam.sp.bed"))
	longsplitreads <- longsplitreads[.distance(longsplitreads, partner(longsplitreads))$max >= minsize * sizemargin]
	longspanningreads <- longspanningreads[.distance(longspanningreads, partner(longspanningreads))$max >= minsize * sizemargin]
}
rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/i/", "~/i/")
#####################################
# Load VCFs
metadata <- LoadMetadata(paste0(rootdir, "data.na12878"))
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.na12878"), metadata=metadata)
vcfs <- sapply(names(vcfs), function(id) {
	gr <- vcfs[[id]]
	seqlevelsStyle(gr) <- "UCSC"
	# only looking at intrachromosomal events at least 50bp in size
	gr <- gr[!is.na(gr$svLen) & abs(gr$svLen) >= minsize,]
	# deletion-like events
	gr <- gr[gr$svtype == "DEL" | (gr$svtype == "BND" & strand(gr) != strand(partner(gr)) & strand(gr) == ifelse(start(gr) < start(partner(gr)), "+", "-")),]
	# on primary chromosomes not overlapping blacklist
	gr <- gr[!overlapsAny(gr, encodeblacklist) & seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")),]
	gr <- gr[gr$partner %in% names(gr),]
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)

#####################################
# Mills
calls_default <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, truthgr=vcfs[["00000000000000000000000000000002"]])
calls_default$calls$CallSet <- "High confidence only"
calls_default$truth$CallSet <- "High confidence only"

calls_all <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, truthgr=vcfs[["00000000000000000000000000000002"]])
calls_all$calls$CallSet <- "High & Low confidence"
calls_all$truth$CallSet <- "High & Low confidence"
mcalls <- rbind(calls_default$calls, calls_all$calls)

# Duplicates are considered FPs
mcalls$fp[mcalls$duptp] <- TRUE
mcalls$tp[mcalls$duptp] <- FALSE

#####################################
# Long read
lrcalls <- rbind_all(lapply(names(vcfs)[names(vcfs) %in% (metadata %>% filter(!is.na(CX_CALLER)))$Id], function(id) {
  write(paste0("Processing ", id), stderr())
  callgr <- vcfs[[id]]
  result <- data.frame(vcfId=callgr$vcfId, QUAL=callgr$QUAL, FILTER=callgr$FILTER) %>% mutate(Id=id, srhits=0, sphits=0)

  .hitCounts <- function(truthgr) {
    hitscounts <- rep(0, length(callgr))
    hits <- findMatchingBreakpoints(callgr, truthgr, maxgap=maxgap, ignore.strand=FALSE, sizemargin=sizemargin)
    hits$QUAL <- callgr$QUAL[hits$queryHits]
    # assign supporting evidence to the call with the highest QUAL
    hits <- hits %>%
      arrange(desc(QUAL), queryHits) %>%
      distinct(subjectHits) %>%
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
  mutate(tp=srhits>=3 | sphits>=7, fp=!tp, CallSet="High & Low confidence")
lrcalls <- rbind(lrcalls, lrcalls %>%
  filter(FILTER %in% c(".", "PASS")) %>%
  mutate(CallSet="High confidence only"))

#####################################
# Plots
.mostSensitivePerCaller <- function(calls, callers=c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra/breakdancer", "cortex", "hydra")) {
  calls %>%
    dplyr::select(Id, CallSet, tp) %>%
    group_by(Id, CallSet) %>%
    summarise(tp=sum(tp)) %>%
    ungroup() %>%
    arrange(desc(tp)) %>%
    left_join(metadata) %>%
    distinct(CallSet, StripCallerVersion(CX_CALLER)) %>%
    filter(StripCallerVersion(CX_CALLER) %in% callers) %>%
    dplyr::select(Id, CallSet)
}
.plotgraphs <- function(calls, label) {
  # use aligner with best sensitivity
  sensAligner <- .mostSensitivePerCaller(calls)
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
    mutate(caller=StripCallerVersion(CX_CALLER), CallSet=relevel(factor(CallSet), "High confidence only"))
  ggplot(roc) +
    aes(group=paste(Id, CallSet), y=tp/2, x=fp/2, linetype=CallSet, color=caller) +
    geom_line() +
    coord_cartesian(xlim=c(0, 1000), ylim=c(0, 3000)) +
    scale_colour_brewer(palette="Set1") +
	scale_linetype_manual(values=c("solid", "dotted")) +
    labs(y="True Positives", x="False Positives")
  saveplot(paste0("na12878_tp_fp_", label, "_error_", maxgap, "bp_", sizemargin, "x"), width=7, height=5)
  ggplot(rbind(roc, roc %>%
  		# add in horizontal line to the y axis
	  	filter(tp > 0) %>%
  		arrange(tp) %>%
	  	distinct(caller, CallSet) %>%
  		mutate(tp=0))) +
    aes(group=paste(Id, CallSet), x=tp/2, y=precision, linetype=CallSet, color=caller) +
    geom_line() +
    coord_cartesian(xlim=c(0, 4000)) +
    scale_colour_brewer(palette="Set1") +
	scale_linetype_manual(values=c("solid", "dotted")) +
    labs(x="True positives", y="Precision")
  saveplot(paste0("na12878_prec_", label, "_error_", maxgap, "bp_", sizemargin, "x"), width=7, height=5)
}

.plotgraphs(mcalls, "Mills")
.plotgraphs(lrcalls, "pacbiomoleculo")

dterror <- mcalls %>%
  dplyr::select(Id, CallSet, bperror) %>%
  left_join(metadata) %>%
  mutate(caller=StripCallerVersion(CX_CALLER))
ggplot(dterror) +
  aes(group=paste(Id, CallSet), x=bperror, color=caller, linetype=CallSet) +
  geom_density() +
	scale_linetype_manual(values=c("solid", "dotted")) +
  scale_colour_brewer(palette="Set1") +
	facet_wrap(~caller) +
  labs(title="Error margin", y="density", x="bp error")




