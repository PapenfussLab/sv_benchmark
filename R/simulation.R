source("sv_benchmark.R")
source("libplot.R")
library(dplyr)
library(stringr)

rootdir <- ifelse(as.character(Sys.info())[1] == "Windows", "W:/i/", "~/i/")

maxgap <- 180
sizemargin <- 0.50
ignore.strand <- TRUE
ignore.duplicates <- TRUE
nominalPosition <- TRUE

#rd <- 60

metadata <- LoadMetadata(paste0(rootdir, "data.rd"))
#metadata <- metadata %>% filter(is.na(CX_CALLER) | (CX_READ_DEPTH==rd))

vcfs <- LoadMinimalSVFromVCF(paste0(rootdir, "data.rd"), metadata=metadata, nominalPosition=nominalPosition)

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


calls_default <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=FALSE, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=TRUE)
mcalls_default <- rbind(calls_default$calls %>% filter(!tp), calls_default$truth %>% mutate(duptp=FALSE))
mcalls_default$CallSet <- "High confidence only"

calls_all <- ScoreVariantsFromTruth(vcfs, metadata, includeFiltered=TRUE, maxgap=maxgap, sizemargin=sizemargin, ignore.strand=TRUE)
mcalls_all <- rbind(calls_all$calls %>% filter(!tp), calls_all$truth %>% mutate(duptp=FALSE))
mcalls_all$CallSet <- "High & Low confidence"
mcalls <- rbind(mcalls_all, mcalls_default)

# Duplicate TP calls
if (ignore.duplicates) {
	# Technically, these events are misclassified
	mcalls <- mcalls %>% filter(!duptp)
}

# Sanity checks
ggplot(mcalls %>%
				 filter(!fp) %>%
				 filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
				 group_by(Id, CallSet, maxgap, ignore.strand, svLen) %>%
				 summarise(sens=sum(tp)/sum(tp+fn)) %>%
				 ungroup() %>%
				 left_join(metadata)) +
		aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, linetype=CX_ALIGNER %na% "", color=CallSet) +
		geom_point() +
		geom_line() +
		scale_x_svlen +
		facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) +
		labs(title="Results sanity check", color="args")
ggplot(mcalls %>%
				 dplyr::select(Id, CallSet, QUAL, tp, fp, fn) %>%
				 filter(Id %in% metadata$Id[str_detect(metadata$CX_REFERENCE_VCF_VARIANTS, "BP")]) %>%
				 arrange(desc(QUAL)) %>%
				 group_by(Id, CallSet) %>%
				 mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
				 group_by(Id, CallSet, QUAL) %>%
				 summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
				 left_join(metadata)) +
		aes(group=paste(Id, CallSet), y=tp, x=fp + 1, linetype=CX_ALIGNER %na% "", color=paste0(CX_CALLER_ARGS %na% "",CX_CALLER_FLAGS %na% "")) +
		geom_line() +
		geom_point(alpha=0.1) +
		facet_grid(CX_CALLER ~ CX_REFERENCE_VCF_VARIANTS) +
		scale_x_log10() +
		labs(title="Results sanity check", color="args")

# use aligner with best sensitivity
sensAligner <- mcalls %>%
	dplyr:::select(Id, CallSet, maxgap, ignore.strand, tp) %>%
	group_by(Id, CallSet, maxgap, ignore.strand) %>%
	summarise(tp=sum(tp)) %>%
	ungroup() %>%
	arrange(desc(tp)) %>%
	left_join(metadata) %>%
	distinct(CallSet, StripCallerVersion(CX_CALLER), CX_READ_LENGTH, CX_READ_DEPTH, CX_READ_FRAGMENT_LENGTH, CX_REFERENCE_VCF, .keep_all = TRUE)

#####################
# Fixed read depth plots
# Separate plot per caller
es <- mcalls %>%
	filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
	filter(!fp) %>%
	filter(Id %in% metadata$Id[metadata$CX_REFERENCE_VCF_VARIANTS %in% c("hetDEL","hetINS","hetINV","hetDUP")]) %>%
	group_by(Id, CallSet, maxgap, ignore.strand, svLen) %>%
	summarise(sens=sum(tp)/sum(tp+fn)) %>%
	ungroup() %>%
	left_join(metadata) %>%
	mutate(caller=StripCallerVersion(CX_CALLER), eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS)) %>%
	filter(caller %in% c("gridss", "breakdancer", "delly", "lumpy", "pindel", "socrates", "manta", "cortex", "hydra"))
roc <- mcalls %>%
	dplyr::select(Id, CallSet, QUAL, tp, fp, fn) %>%
	rbind(mcalls %>% dplyr::select(Id, CallSet) %>% distinct(Id, CallSet) %>% mutate(QUAL=max(mcalls$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
	filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
	group_by(Id, CallSet) %>%
	arrange(desc(QUAL)) %>%
	mutate(events=sum(tp) + sum(fn)) %>%
	mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
	group_by(Id, CallSet, QUAL, events) %>%
	summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
	ungroup() %>%
	mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events) %>%
	left_join(metadata) %>%
	mutate(caller=StripCallerVersion(CX_CALLER), eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS)) %>%
	filter(caller %in% c("gridss", "breakdancer", "delly", "lumpy", "pindel", "socrates", "manta", "cortex", "hydra"))
# reduction of roc plot points by elimination of points on straight lines
roc <- roc %>%
	group_by(Id, CallSet) %>%
	arrange(tp + fp) %>%
	filter(
		# keep start/end
		is.na(lag(tp)) | is.na(lead(tp)) |
			# keep group transitions (TODO: is there a way to make lead/lag across group_by return NA?)
			Id != lag(Id) | CallSet != lag(CallSet) |
			Id != lead(Id) | CallSet != lead(CallSet) |
			# slopes not equal dx1/dy1 != dx2/dy2 -> dx1*dy2 != dx2*dy1
			(tp - lag(tp))*(lead(fp) - lag(fp)) != (lead(tp) - lag(tp))*(fp - lag(fp)) |
			# less than 10 calls wide
			lead(tp) - lag(tp) > 10 |
			# keep every 5th row
			row_number() %% 10 == 0) %>%
	ungroup()
####
# rd
ggplot(es) +
	aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, color=eventtype, linetype=CallSet) +
	geom_line(size=0.5) +
	scale_x_svlen +
	facet_grid(caller ~ CX_READ_DEPTH, switch="y") +
	labs(title="", y="Sensitivity", x="Event size", color="Event Type", linetype="Call set")
saveplot(paste0("sim_per_caller_event_size_line"), width=600, height=300, units=c("mm"))
write.csv(es, "sim_per_caller_event_size_line.csv")
ggplot(roc %>% arrange(desc(QUAL))) +
	aes(group=paste(Id, CallSet), y=sens, x=fp+1, linetype=CallSet, color=eventtype) +
	geom_line() +
	scale_x_continuous(breaks=c(1, 11, 101, 1001, 10001, 100001),
										 labels=c("0", "10", "100", "1k", "10k", "100k"),
										 minor_breaks=c(),
										 trans="log10") +
	facet_grid(caller ~ CX_READ_DEPTH) +
	labs(title="", y="Sensitivity", x="False Positives", color="Call set", linetype="Call set")
saveplot(paste0("sim_per_caller_roc"), width=300, height=300, units=c("mm"))
write.csv(roc %>% arrange(desc(QUAL)), "sim_per_caller_roc.csv")
for (rd in unique(metadata$CX_READ_DEPTH)) {
	ggplot(es %>% filter(CX_READ_DEPTH==rd)) +
		aes(group=paste(Id, CallSet), x=abs(svLen), y=sens) +
		geom_area(aes(fill=CallSet), position="identity") +
		geom_area(data=es[es$CallSet=="High confidence only",], aes(fill=CallSet), position="identity") +
		scale_fill_manual(values=c("#2166ac", "#67a9cf")) +
		scale_x_svlen_short +
		facet_grid(caller ~ eventtype, switch="y") +
		labs(title="", y="Sensitivity", x="Event size", color="Call set", linetype="Call set")
	saveplot(paste0("sim_", rd, "x_per_caller_event_size_fill"), width=220, height=300, units=c("mm"))
	saveplot(paste0("ppt_sim_", rd, "x_per_caller_event_size_fill"), width=220*16/10, height=220, units=c("mm"))
	write.csv(es, paste0("sim_", rd, "x_per_caller_event_size.csv"))
	ggplot(es %>% filter(CX_READ_DEPTH==rd)) +
		aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, color=eventtype, linetype=CallSet) +
		geom_line(size=0.5) +
		scale_x_svlen +
		facet_grid(caller ~ ., switch="y") +
		labs(title="", y="Sensitivity", x="Event size", color="Event Type", linetype="Call set")
	saveplot(paste0("sim_", rd, "x_per_caller_event_size_line"))
	plotdata <- roc %>% filter(CX_READ_DEPTH==rd) %>% arrange(desc(QUAL))
	ggplot(plotdata) +
		aes(group=paste(Id, CallSet), y=sens, x=fp+1, color=CallSet) +
		scale_color_manual(values=c("#2166ac", "#67a9cf")) +
		#scale_color_brewer(palette="Paired") +
		geom_line(size=0.2) +
		geom_point(size=0.5) +
		geom_point(size=0.5, data=plotdata %>% filter(CallSet=="High confidence only")) +
		scale_x_continuous(breaks=c(1, 11, 101, 1001, 10001),
											 labels=c("", "10", "100", "1k", "10k"),
											 minor_breaks=c(),
											 trans="log10") +
		facet_grid(caller ~ eventtype, switch="y") +
		labs(title="", y="Sensitivity", x="False Positives", color="Call set", linetype="Call set")
	saveplot(paste0("sim_", rd, "x_per_caller_roc"), width=300, height=300, units=c("mm"))
	saveplot(paste0("ppt_sim_", rd, "x_per_caller_roc"), width=220*16/10, height=220, units=c("mm"))
	write.csv(roc %>% arrange(desc(QUAL)), paste0("sim_", rd, "x_per_caller_roc.csv"))
}
# table of maximum sensitivity
roc %>%
	dplyr::select(caller, CallSet, eventtype, sens, fp) %>%
	group_by(caller, CallSet, eventtype) %>%
	summarise(sens=max(sens), fp=max(fp)) %>%
	ungroup() %>%
	arrange(eventtype, desc(sens))

# table of event size detection range (Supp table 1b)
write.csv(es %>%
  group_by(caller, eventtype, svLen, sens) %>%
  group_by(caller, eventtype, svLen) %>%
  summarise(sens=max(sens)) %>%
  ungroup() %>%
  mutate(svLen=abs(svLen)) %>%
  arrange(caller, eventtype, svLen),
  "event_size_sens.csv")

######
# Supp Fig
######
# positional delta between called position and truth position
errSummary <- mcalls %>%
  filter(abs(svLen) > 50 | is.na(svLen)) %>%
  filter(tp) %>%
  filter(CallSet=="High & Low confidence") %>%
  filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
  left_join(metadata) %>%
  mutate(caller=StripCallerVersion(CX_CALLER), eventtype=PrettyVariants(CX_REFERENCE_VCF_VARIANTS)) %>%
  filter(caller %in% c("gridss", "breakdancer", "delly", "lumpy", "pindel", "socrates", "manta", "cortex", "hydra"))

ggplot() +
  aes(x=error, fill=positionLabel) +
  facet_grid(eventtype ~ caller) +
  geom_histogram(data=errSummary %>% mutate(error=bperror, positionLabel="Nominal position"), binwidth=2, position="identity", na.rm=TRUE) +
  geom_histogram(data=errSummary %>% mutate(error=pmax(0, bperror-HOMLEN), positionLabel="Allowing for\nreported micro-homology"), binwidth=2, position="identity", na.rm=TRUE) +
  scale_y_continuous(breaks=10**(0:4), labels=10**(0:4), trans="log1p") +
  #scale_y_log10() +
  scale_fill_brewer(type="qual", palette="Dark2", name="Error margin") + 
  labs(x="Error (bp)", y="Events") +
  theme_bw()
saveplot(paste0("sim_bp_error"), width=300, height=220, units=c("mm"))

# table of maximum F1-score (Supp table 1b)
roc %>%
  dplyr::mutate(f1score=2*precision*sens/(precision+sens)) %>%
  dplyr::select(caller, eventtype, f1score) %>%
  dplyr::filter(!is.na(f1score)) %>%
  group_by(caller, eventtype) %>%
  summarise(f1score=max(f1score)) %>%
  ungroup() %>%
  arrange(eventtype, caller) %>%
  filter(eventtype %in% c("SINE/ALU Breakpoint", "Breakpoint"))

roccombined <- mcalls %>%
	dplyr::select(Id, CallSet, QUAL, tp, fp, fn) %>%
	rbind(mcalls %>% dplyr::select(Id, CallSet) %>% distinct(Id, CallSet) %>% mutate(QUAL=max(mcalls$QUAL) + 1, tp=0, fp=0, fn=0)) %>%
	filter(paste(Id, CallSet) %in% paste(sensAligner$Id, sensAligner$CallSet)) %>%
	left_join(metadata %>% mutate(caller=StripCallerVersion(CX_CALLER)) %>% dplyr::select(Id, caller)) %>%
	filter(caller %in% c("gridss", "breakdancer", "delly", "lumpy", "pindel", "socrates", "manta", "cortex", "hydra")) %>%
	group_by(caller, CallSet) %>%
	arrange(desc(QUAL)) %>%
	mutate(events=sum(tp) + sum(fn)) %>%
	mutate(tp=cumsum(tp), fp=cumsum(fp), fn=cumsum(fn)) %>%
	group_by(caller, CallSet, QUAL, events) %>%
	summarise(tp=max(tp), fp=max(fp), fn=max(fn)) %>%
	ungroup() %>%
	mutate(precision=tp / (tp + fp), fdr=1-precision, sens=tp/events)
ggplot(roccombined %>% arrange(desc(QUAL))) +
	aes(group=paste(caller, CallSet), y=sens, x=fp+1, linetype=CallSet, color=CallSet) +
	scale_color_manual(values=c("#000000", "#666666")) +
	geom_line() +
	scale_x_continuous(breaks=c(1, 11, 101, 1001, 10001, 100001),
										 labels=c("", "10", "100", "1k", "10k", "100k"),
										 minor_breaks=c(),
										 trans="log10") +
	facet_wrap(~ caller) +
	labs(title="", y="Sensitivity", x="False Positives", color="Call set", linetype="Call set")
saveplot(paste0("sim_", rd, "x_per_caller_roc_combined"), width=300, height=300, units=c("mm"))
write.csv(roccombined %>% arrange(desc(QUAL)), paste0("sim_", rd, "x_per_caller_roc_combined.csv"))
