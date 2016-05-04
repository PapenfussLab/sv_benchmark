source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)

theme_set(theme_bw())

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
	distinct(Filter, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF)

sens <- mcalls %>%
	filter(!fp) %>%
	filter(Id %in% sensAligner$Id) %>%
	filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	group_by(Id, Filter, maxgap, ignore.strand, Filter, svLen) %>%
	summarise(sens=sum(tp)/sum(tp+fn)) %>%
	ungroup() %>%
	left_join(metadata) %>%
	mutate(caller=StripCallerVersion(CX_CALLER)) %>%
	filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra", "cortex", "hydra"))

for (rl in unique(sens$CX_READ_LENGTH)) {
for (rd in unique(sens$CX_READ_DEPTH)) {
for (fragsize in unique(sens$CX_READ_FRAGMENT_LENGTH)) {
	ggplot(sens %>% filter(CX_READ_DEPTH==rd & CX_READ_FRAGMENT_LENGTH==fragsize & CX_READ_LENGTH==rl)) +
		aes(group=Id, x=abs(svLen), y=sens^5, group=Id, color=CX_CALLER, shape=caller, linetype=CX_ALIGNER) +
		geom_point() +
		geom_line() + scale_x_log10() +
		scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")) +
		facet_grid(Filter ~ PrettyVariants(CX_REFERENCE_VCF_VARIANTS))
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
		filter(str_detect(CX_REFERENCE_VCF_VARIANTS, "BP")) %>%
		mutate(caller=StripCallerVersion(CX_CALLER)) %>%
		filter(caller %in% c("delly"))
	) +
    aes(y=tp, x=fp, color=CX_CALLER, linetype=CX_ALIGNER) +
    geom_line() +
    geom_point() +
    facet_grid(Filter ~ CX_REFERENCE_VCF_VARIANTS)



