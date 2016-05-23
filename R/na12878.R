source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)

#rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/i/", "~/i/")

maxgap <- 180
sizemargin <- 0.5
ignore.strand <- TRUE

#####################################
# hg19 blacklist
encodeblacklist <- import(paste0(rootdir, "input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed"))

#####################################
# na12878 Moleculo/PacBio truth
if (!exists("srpacbio") || !exists("srmoleculo") || !exists("sppacbio") || !exists("spmoleculo")) {
	.loadlongreadbed <- function(filename) {
		lrbed <- import.bed(con=paste0(rootdir, filename))
		seqlevelsStyle(lrbed) <- "UCSC"
		return(lrbed)
	}
	srpacbio <- c(
		.loadlongreadbed("input.na12878/chemistry1.sorted.bam.sr.bam.sr.bed"),
		.loadlongreadbed("input.na12878/chemistry_2_picard.bam.sr.bam.sr.bed"),
		.loadlongreadbed("input.na12878/chemistry_3_picard.bam.sr.bam.sr.bed"),
		.loadlongreadbed("input.na12878/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sr.bam.sr.bed"))
	srmoleculo <- .loadlongreadbed("input.na12878/NA12878.moleculo.bwa-mem.20140110.bam.sr.bam.sr.bed")
	sppacbio <- c(
		.loadlongreadbed("input.na12878/chemistry1.sorted.bam.sp.bed"),
		.loadlongreadbed("input.na12878/chemistry_2_picard.bam.sp.bed"),
		.loadlongreadbed("input.na12878/chemistry_3_picard.bam.sp.bed"),
		.loadlongreadbed("input.na12878/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sp.bed"))
	spmoleculo <- .loadlongreadbed("input.na12878/NA12878.moleculo.bwa-mem.20140110.bam.sp.bed")
}
#####################################
# Load VCFs
metadata <- LoadMetadata(paste0(rootdir, "data.na12878"))
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.na12878"), metadata=metadata)
vcfs <- sapply(names(vcfs), function(id) {
	gr <- vcfs[[id]]
	seqlevelsStyle(gr) <- "UCSC"
	# only looking at intrachromosomal events at least 50bp in size
	gr <- gr[!is.na(gr$svLen) & abs(gr$svLen) > 50,]
	# deletion-like events
	gr <- gr[gr$svtype == "DEL" | (gr$svtype == "BND" & strand(gr) != strand(partner(gr)) & strand(gr) == ifelse(start(gr) < start(partner(gr)), "+", "-")),]
	# on primary chromosomes not overlapping blacklist
	gr <- gr[!overlapsAny(gr, encodeblacklist) & seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")),]
	gr <- gr[gr$partner %in% names(gr),]
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)



calls_default <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, truthgr=vcfs[["00000000000000000000000000000002"]])
mcalls_default <- rbind(calls_default$calls %>% filter(!tp), calls_default$truth)
mcalls_default$CallSet <- "High & Low confidence"

calls_all <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, truthgr=vcfs[["00000000000000000000000000000002"]])
mcalls_all <- rbind(calls_all$calls %>% filter(!tp), calls_all$truth)
mcalls_all$CallSet <- "High confidence only"
mcalls <- rbind(mcalls_default, mcalls_all)

# use aligner with best sensitivity
sensAligner <- mcalls %>%
	select(Id, CallSet, maxgap, ignore.strand, tp) %>%
	group_by(Id, CallSet, maxgap, ignore.strand) %>%
	summarise(tp=sum(tp)) %>%
	ungroup() %>%
	arrange(desc(tp)) %>%
	left_join(metadata) %>%
	distinct(CallSet, StripCallerVersion(CX_CALLER), CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH)

roc <- mcalls %>%
    filter(Id %in% sensAligner$Id) %>%
    select(Id, CallSet, QUAL, tp, fp, fn) %>%
    arrange(desc(QUAL)) %>%
    group_by(Id, CallSet) %>%
    mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
    group_by(Id, CallSet, QUAL) %>%
    summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
    mutate(precision=tp / (tp + fp), fdr=1-precision) %>%
    left_join(metadata) %>%
    mutate(caller=StripCallerVersion(CX_CALLER))
ggplot(roc )+#%>% filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra/breakdancer", "cortex", "hydra"))) +
  aes(group=interaction(Id, CallSet), y=tp/2, x=fp/2, linetype=CallSet, color=CX_CALLER) +
  geom_line() +
  scale_x_continuous(limits=c(0, 1000)) +
  labs(title="ROC per caller")



