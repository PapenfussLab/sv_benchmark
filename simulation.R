source("sv_benchmark.R")
library(ggplot2)

maxgap <- 200
ignore.strand <- TRUE

metadata <- LoadMetadata("data.rd")
vcfs <- LoadMinimalSVFromVCF("data.rd", metadata=metadata)

calls_default <- ScoreVariantsFromTruthVCF(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, ignore.strand=TRUE)
mcalls_default <- rbind(calls_default$calls %>% filter(!tp), calls_default$truth)
mcalls_default$filter <- "Default calls"
calls_all <- ScoreVariantsFromTruthVCF(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, ignore.strand=TRUE)
mcalls_all <- rbind(calls_all$calls %>% filter(!tp), calls_all$truth)
mcalls_all$filter <- "All calls"
mcalls <- rbind(mcalls_default, mcalls_all)

# use aligner with best sensitivity
sensAligner <- mcalls %>%
	select(Id, filter, tp) %>%
	group_by(Id, filter) %>%
	summarise(tp=sum(tp)) %>%
	arrange(desc(tp)) %>%
	left_join(metadata) %>%
	distinct(filter, CX_CALLER, CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF)

sens <- mcalls %>%
	filter(!fp) %>%
	filter(Id %in% sensAligner$Id) %>%
	filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	group_by(Id, svLen) %>%
	summarise(sens=sum(tp)/sum(tp+fn)) %>%
	left_join(metadata) %>%
	mutate(caller=str_extract(CX_CALLER, "^[^/]+")) %>%
	filter(caller %in% c("gridss", "breakdancer", "cortex", "delly", "lumpy", "pindel", "socrates"))

ggplot(sens %>% filter(CX_READ_DEPTH==60)) +
	aes(x=abs(svLen), y=sens^5, group=Id, color=CX_CALLER, shape=caller, linetype=caller) +
	geom_point() +
	geom_line() + scale_x_log10() +
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)^5, labels=c("0", "", "", "", "", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")) +
	facet_wrap(~CX_REFERENCE_VCF_VARIANTS)
