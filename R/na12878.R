source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)

#rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/i/", "~/i/")

maxgap <- 200
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
	seqlevelsStyle(gr) <- "UCSC"
	gr <- vcfs[[id]]
	# only looking at intrachromosomal events at least 50bp in size
	gr <- gr[!is.na(gr$svLen) & abs(gr$svLen) > 50,]
	# deletion-like events
	gr <- gr[gr$svtype == "DEL" | (gr$svtype == "BND" & strand(gr) != strand(partner(gr)) & strand(gr) == ifelse(start(gr) < start(partner(gr)), "+", "-")),]
	# on primary chromosomes not overlapping blacklist
	gr <- gr[!overlapsAny(gr, encodeblacklist) & seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")),]
	gr <- gr[gr$partner %in% names(gr),]
	return(gr)
}, simplify=FALSE, USE.NAMES=TRUE)



#calls_default <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=TRUE, truthgr=vcfs[["00000000000000000000000000000002"]])
#mcalls_default <- rbind(calls_default$calls %>% filter(!tp), calls_default$truth)
#mcalls_default$Filter <- "Default calls"
mcalls_default <- NULL

calls_all <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=TRUE, truthgr=vcfs[["00000000000000000000000000000002"]])
mcalls_all <- rbind(calls_all$calls %>% filter(!tp), calls_all$truth)
mcalls_all$Filter <- "All calls"
mcalls <- rbind(mcalls_default, mcalls_all)

# use aligner with best sensitivity
sensAligner <- mcalls %>%
	select(Id, Filter, maxgap, ignore.strand, tp) %>%
	group_by(Id, Filter, maxgap, ignore.strand) %>%
	summarise(tp=sum(tp)) %>%
	ungroup() %>%
	arrange(desc(tp)) %>%
	left_join(metadata) %>%
	distinct(Filter, StripCallerVersion(CX_CALLER), CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH)

ggplot(mcalls %>%
    filter(Id %in% sensAligner$Id) %>%
    select(Id, Filter, QUAL, tp, fp, fn) %>%
    arrange(desc(QUAL)) %>%
    group_by(Id, Filter) %>%
    mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
    group_by(Id, Filter, QUAL) %>%
    summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
    mutate(precision=tp / (tp + fp), fdr=1-precision) %>%
    left_join(metadata) %>%
    mutate(caller=StripCallerVersion(CX_CALLER)) %>%
    filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra/delly", "cortex", "hydra"))) +
  aes(group=Id, y=tp, x=fp, linetype=Filter, color=caller) +
  geom_line() +
  geom_point() +
  labs(title="ROC per caller")



# Comparison of GRIDSS parameters
# mcalls %>%
#   select(Id, Filter, tp) %>%
#   group_by(Id, Filter) %>%
#   summarise(tp=sum(tp)) %>%
#   ungroup() %>%
#   left_join(metadata) %>%
#   select(tp, Id, Filter, CX_REFERENCE_VCF_VARIANTS, CX_CALLER, CX_CALLER_ARGS) %>%
#   filter(str_detect(CX_CALLER, "gridss")) %>%
#   arrange(CX_REFERENCE_VCF_VARIANTS, desc(tp)) %>%
#   filter(Filter=="All calls") %>%
#   View

sens <- mcalls %>%
	filter(!fp) %>%
	filter(Id %in% sensAligner$Id) %>%
	filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	group_by(Id, Filter, maxgap, ignore.strand, svLen) %>%
	summarise(sens=sum(tp)/sum(tp+fn)) %>%
	ungroup() %>%
	left_join(metadata) %>%
	mutate(caller=StripCallerVersion(CX_CALLER)) %>%
  filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra/delly", "cortex", "hydra"))
	#filter(caller %in% c("gridss", "breakdancer", "tigra/breakdancer", "cortex", "hydra"))

for (rl in unique(sens$CX_READ_LENGTH)) {
for (rd in unique(sens$CX_READ_DEPTH)) {
for (fragsize in unique(sens$CX_READ_FRAGMENT_LENGTH)) {
	ggplot(sens %>% filter(CX_READ_DEPTH==rd & CX_READ_FRAGMENT_LENGTH==fragsize & CX_READ_LENGTH==rl)) +
		aes(group=Id, x=abs(svLen), y=sens^5, color=CX_CALLER, shape=caller, linetype=CX_ALIGNER) +
		geom_point() +
		geom_line() +
    scale_x_svlen +
		scale_y_power5 +
		facet_grid(Filter ~ CX_REFERENCE_VCF_VARIANTS)
}}}

roc <- mcalls %>%
	select(Id, Filter, QUAL, tp, fp, fn) %>%
	arrange(desc(QUAL)) %>%
	group_by(Id, Filter) %>%
	mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
	group_by(Id, Filter, QUAL) %>%
	summarise(tp=max(tp), fp=max(fp), fn=max(fn))

for (var in unique(metadata$CX_REFERENCE_VCF_VARIANTS)) {
  ggplot(roc %>%
  		left_join(metadata) %>%
  		mutate(caller=StripCallerVersion(CX_CALLER)) %>%
  		filter(Filter=="All calls") %>%
		  filter(CX_REFERENCE_VCF_VARIANTS==var)
  	) +
      aes(group=Id, y=tp, x=fp + 1, color=CX_CALLER, linetype=CX_ALIGNER, shape=Filter) +
      geom_line() +
      geom_point(alpha=0.1) +
      facet_wrap(~ CX_CALLER) +
      scale_x_log10() +
      labs(title=var)
}


