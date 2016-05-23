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


calls_default <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=TRUE)
mcalls_default <- rbind(calls_default$calls %>% filter(!tp), calls_default$truth)
mcalls_default$CallSet <- "High & Low confidence"
mcalls_default <- NULL

calls_all <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=TRUE)
mcalls_all <- rbind(calls_all$calls %>% filter(!tp), calls_all$truth)
mcalls_all$CallSet <- "High confidence only"
mcalls <- rbind(mcalls_default, mcalls_all)

# Sanity checks
ggplot(mcalls %>%
         filter(!fp) %>%
         filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
         group_by(Id, CallSet, maxgap, ignore.strand, svLen) %>%
         summarise(sens=sum(tp)/sum(tp+fn)) %>%
         ungroup() %>%
         left_join(metadata)) +
    aes(group=Id, x=abs(svLen), y=sens, linetype=CX_ALIGNER %na% "", color=paste0(CX_CALLER_ARGS %na% "",CX_CALLER_FLAGS %na% "")) +
    geom_point() +
    geom_line() +
    scale_x_svlen +
    facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) + 
    labs(title="Results sanity check", color="args")

ggplot(mcalls %>%
         select(Id, CallSet, QUAL, tp, fp, fn) %>%
         filter(Id %in% metadata$Id[str_detect(metadata$CX_REFERENCE_VCF_VARIANTS, "BP")]) %>%
         arrange(desc(QUAL)) %>%
         group_by(Id, CallSet) %>%
         mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
         group_by(Id, CallSet, QUAL) %>%
         summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
         left_join(metadata)) +
    aes(group=Id, y=tp, x=fp + 1, linetype=CX_ALIGNER %na% "", color=paste0(CX_CALLER_ARGS %na% "",CX_CALLER_FLAGS %na% "")) +
    geom_line() +
    geom_point(alpha=0.1) +
    facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) +
    scale_x_log10() + 
    labs(title="Results sanity check", color="args")

# use aligner with best sensitivity
sensAligner <- mcalls %>%
	select(Id, CallSet, maxgap, ignore.strand, tp) %>%
	group_by(Id, CallSet, maxgap, ignore.strand) %>%
	summarise(tp=sum(tp)) %>%
	ungroup() %>%
	arrange(desc(tp)) %>%
	left_join(metadata) %>%
	distinct(CallSet, StripCallerVersion(CX_CALLER), CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF)

# Separate plot per caller
ggplot(mcalls %>%
         filter(Id %in% sensAligner$Id) %>%
         filter(!fp) %>%
         filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
         group_by(Id, CallSet, maxgap, ignore.strand, svLen) %>%
         summarise(sens=sum(tp)/sum(tp+fn)) %>%
         ungroup() %>%
         left_join(metadata) %>%
         mutate(caller=StripCallerVersion(CX_CALLER)) %>%
         filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra/delly", "cortex", "hydra"))) +
  aes(group=Id, x=abs(svLen), y=sens, shape=CX_REFERENCE_VCF_VARIANTS, color=CX_REFERENCE_VCF_VARIANTS, linetype=CallSet) +
  scale_shape_manual(values=c(2,3,4,6)) + 
  geom_point(size=2) +
  geom_line(size=0.5, alpha=0.5) +
  scale_x_svlen +
  facet_grid(caller ~ .) + 
  labs(title="sensitivity per caller")

ggplot(mcalls %>%
    filter(Id %in% sensAligner$Id) %>%
    select(Id, CallSet, QUAL, tp, fp, fn) %>%
    arrange(desc(QUAL)) %>%
    group_by(Id, CallSet) %>%
    mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn), events=sum(tp) + sum(fn)) %>%
    group_by(Id, CallSet, QUAL, events) %>%
    summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
    mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
    left_join(metadata) %>%
    mutate(caller=StripCallerVersion(CX_CALLER)) %>%
    filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates", "tigra/delly", "cortex", "hydra"))) + 
  aes(group=Id, y=fdr, x=sens, linetype=CallSet) +
  geom_line() +
  facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) +
  labs(title="ROC per caller")



# Comparison of GRIDSS parameters
# mcalls %>%
#   select(Id, CallSet, tp) %>%
#   group_by(Id, CallSet) %>%
#   summarise(tp=sum(tp)) %>%
#   ungroup() %>%
#   left_join(metadata) %>%
#   select(tp, Id, CallSet, CX_REFERENCE_VCF_VARIANTS, CX_CALLER, CX_CALLER_ARGS) %>%
#   filter(str_detect(CX_CALLER, "gridss")) %>%
#   arrange(CX_REFERENCE_VCF_VARIANTS, desc(tp)) %>%
#   filter(CallSet=="All calls") %>%
#   View

sens <- mcalls %>%
	filter(!fp) %>%
	filter(Id %in% sensAligner$Id) %>%
	filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	group_by(Id, CallSet, maxgap, ignore.strand, svLen) %>%
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
		facet_grid(CallSet ~ CX_REFERENCE_VCF_VARIANTS)
}}}

roc <- mcalls %>%
	select(Id, CallSet, QUAL, tp, fp, fn) %>%
	arrange(desc(QUAL)) %>%
	group_by(Id, CallSet) %>%
	mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
	group_by(Id, CallSet, QUAL) %>%
	summarise(tp=max(tp), fp=max(fp), fn=max(fn))

for (var in unique(metadata$CX_REFERENCE_VCF_VARIANTS)) {
  ggplot(roc %>%
  		left_join(metadata) %>%
  		mutate(caller=StripCallerVersion(CX_CALLER)) %>%
  		filter(CallSet=="All calls") %>%
		  filter(CX_REFERENCE_VCF_VARIANTS==var)
  	) +
      aes(group=Id, y=tp, x=fp + 1, color=CX_CALLER, linetype=CX_ALIGNER, shape=CallSet) +
      geom_line() +
      geom_point(alpha=0.1) +
      facet_wrap(~ CX_CALLER) +
      scale_x_log10() + 
      labs(title=var)
}


