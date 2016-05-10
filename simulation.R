source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/projects/sv_benchmark/", "~/projects/sv_benchmark/")

maxgap <- 200
ignore.strand <- TRUE

metadata <- LoadMetadata(paste0(rootdir, "data.rd"))
metadata <- metadata %>% filter(is.na(CX_CALLER) | (CX_READ_DEPTH==60))
vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.rd"), metadata=metadata)

# matching transforms made in chromothripsis_sim.R
vcfs <- sapply(names(vcfs), function(id) {
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
}, simplify=FALSE, USE.NAMES=TRUE)

calls_default <- ScoreVariantsFromTruthVCF(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=TRUE)
mcalls_default <- rbind(calls_default$calls %>% filter(!tp), calls_default$truth)
mcalls_default$Filter <- "Default calls"
calls_all <- ScoreVariantsFromTruthVCF(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=TRUE)
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
	distinct(Filter, StripCallerVersion(CX_CALLER), CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF)

sens <- mcalls %>%
	filter(!fp) %>%
	filter(Id %in% sensAligner$Id) %>%
	filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	group_by(Id, Filter, maxgap, ignore.strand, Filter, svLen) %>%
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

ggplot(roc %>%
		left_join(metadata) %>%
		mutate(caller=StripCallerVersion(CX_CALLER)) %>%
		filter(Filter=="All calls")
	) +
    aes(group=Id, y=tp, x=fp + 1, color=CX_CALLER, linetype=CX_ALIGNER, shape=Filter) +
    geom_line() +
    geom_point(size=0.2) +
    facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) +
    scale_x_log10()


