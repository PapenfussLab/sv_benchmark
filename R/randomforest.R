library(randomForest)
source("libmanuscript_figures.R")

ignore.interchromosomal <- TRUE
mineventsize <- 51
maxeventsize <- NULL #TODO: what max event size should we use? Need one for chm since truth only has calls under a certain size
maxgap <- 200
sizemargin <- 0.25
ignore.strand <- TRUE
nominalPosition <- FALSE
longread_requiredHits <-5
longread_minMapq <- 0


datadir <- "../data.chm13"
sample_name <- "chm13"
ids <- c(
    "432ca0d47c12f8a3b9d1eb04615ff737" # gridss/1.6.1-SNAPSHOT
    # "c639766990fcfbca2eb45f7806362fe6" # bcftools/1.3.1
	)
truth_id <- "00000000000000000000000000000013"
truth_name <- "Huddleston et al"
grtransformName <- "None"
eventtype <- "DEL"
setCacheRootPath(datadir)
fileprefix <- str_replace(paste(sample_name, truth_name, ifelse(use_roc_fdr, "fdr", ""), paste0(eventtype, collapse = "_"), sep="_"), "[ /]", "_")
all_ids <- c(truth_id, ids)
metadata <- LoadCachedMetadata(datadir)
metadata <- metadata %>% filter(Id %in% all_ids)
# force truth
metadata$CX_REFERENCE_VCF <- list.files(datadir, pattern=paste0("^", truth_id, ".*.vcf$"))
callgr_chm13 <- .CachedLoadCallMatrixForIds(
		datadir=datadir,
		metadata=metadata,
		ids=all_ids,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand,
		grtransform=lroptions$grtransform[[grtransformName]],
		grtransformName=grtransformName,
		nominalPosition=nominalPosition,
		eventtype=eventtype)

filter_to_variables <- function(df) {
	df %>%
		mutate(tp=factor(tp, levels=c("fp", "tp"))) %>%
		mutate(repeatClass=as.factor(repeatClass)) %>%
		select(tp, QUAL, svtype, svLen, insLen, HOMLEN, repeatClass, snp50bp)
}
traindf <- as.data.frame(callgr_chm13) %>%
	filter(!ignore.filtered) %>%
	filter(Id != "00000000000000000000000000000013") %>%
	mutate(tp=ifelse(fId00000000000000000000000000000013 >=0, "tp", "fp")) %>%
	filter_to_variables()

predictors <- traindf %>% select(-tp, -svtype)
chm13_rf <- randomForest(x = predictors, y=traindf$tp, ntree = 2000, nodesize = 5, maxnodes = 16)






datadir <- "../data.chm1"
sample_name <- "chm1"
ids <- c(
    "b1112f1c3cbd28c464f58fc5c5c02f9b" # gridss/1.4.1
    )  # socrates/1.13.1
truth_id <- "00000000000000000000000000000001"
truth_name <- "Huddleston et al"
grtransformName <- "None"
setCacheRootPath(datadir)
fileprefix <- str_replace(paste(sample_name, truth_name, ifelse(use_roc_fdr, "fdr", ""), paste0(eventtype, collapse = "_"), sep="_"), "[ /]", "_")
all_ids <- c(truth_id, ids)
metadata <- LoadCachedMetadata(datadir)
metadata <- metadata %>% filter(Id %in% all_ids)
callgr_chm1 <- .CachedLoadCallMatrixForIds(
		datadir=datadir,
		metadata=metadata,
		ids=all_ids,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand,
		grtransform=lroptions$grtransform[[grtransformName]],
		grtransformName=grtransformName,
		nominalPosition=nominalPosition,
		eventtype=eventtype)
callgr_chm1$truthQUAL <- mcols(callgr_chm1)[,paste0("fId",truth_id)]

testdf_chm1 <- as.data.frame(callgr_chm1) %>%
	filter(!ignore.filtered) %>%
	filter(Id != "00000000000000000000000000000001") %>%
	mutate(tp=ifelse(fId00000000000000000000000000000001 >=0, "tp", "fp")) %>%
	filter_to_variables()

chm1_predictors <- testdf_chm1 %>% select(-tp, -svtype)
mprob_chm1 <- predict(chm13_rf, chm1_predictors, type="prob")

rocdf_chm1 <- as.data.frame(callgr_chm1) %>%
	filter(ignore.filtered) %>%
	filter(Id != truth_id) %>%
	mutate(tp=ifelse(fId00000000000000000000000000000001 >=0, "tp", "fp")) %>%
	filter_to_variables() %>%
	mutate(scoring="gridss") %>%
	select(tp, scoring, QUAL) %>%
	rbind(testdf_chm1 %>% mutate(scoring="gridss (all)") %>%
		select(tp, scoring, QUAL)
	) %>%
	rbind(testdf_chm1 %>%
			select(tp) %>%
			mutate(scoring="random forest", QUAL=mprob_chm1[,2])
	) %>%
	mutate(tp=tp == "tp") %>%
	mutate(fp=!tp) %>%
	group_by(scoring, QUAL) %>%
	summarise(fp=sum(fp), tp=sum(tp)) %>%
	group_by(scoring) %>%
	arrange(desc(QUAL)) %>%
	mutate(
		fp=cumsum(fp),
		tp=cumsum(tp)) %>%
	# each QUAL score is a point on the ROC plott
	group_by(scoring, QUAL) %>%
	summarise(tp=max(tp), fp=max(fp)) %>%
	# QUAL scores with the same number of tp calls can be merged on the ROC plot
	group_by(scoring, tp) %>%
	summarise(fp=max(fp), QUAL=min(QUAL)) %>%
	# subsample along tp and tp+fp axis
	group_by(scoring) %>%
	dplyr::slice(unique(c(
		1,
		findInterval(seq(0, max(tp), max(tp)/1000), tp),
		findInterval(seq(0, max(tp + fp), max(tp + fp)/1000), tp + fp),
		n()
	))) %>%
	mutate(is_endpoint = tp == max(tp)) %>%
	ungroup() %>%
	#TODO: fn and sens need eventCount using group_by(Id, CallSet, !!!groupingCols)
	mutate(
		precision=tp / (tp + fp),
		fdr=1-precision)
randomforest_plot <- ggplot(rocdf_chm1) +
	aes(x=tp, y=precision, colour=scoring) +
	geom_line() +
	geom_point(data = rocdf_chm1 %>% filter(is_endpoint), size = 2) +
	caller_colour_scheme +
	coord_cartesian(ylim = c(0,1)) +
	scale_y_continuous(labels = scales::percent) +
	theme_cowplot() +
	background_grid("xy", "none") +
	labs(
		color = "caller",
		linetype = "call set",
		x = "Number of true positives called",
		y="Precision (1-FDR)",
		title="CHM1 trained on CHM13 sequence context"
		)
saveplot("gridsss_random_forest", plot=randomforest_plot, height=6, width=7)

# NA12878
datadir <- "../data.na12878"
sample_name <- "NA12878"
ids <- c(
	"8aaf2886ffe782e666661d6b890b330a" #gridss/1.4.1
)
grtransformName <- "DAC"
longreadbedpedir <- paste0(datadir, "/../", "input.na12878/longread/")
truth_id <- "00000000000000000000000000000003"
truth_name <- "Parikh et al"
setCacheRootPath(datadir)
fileprefix <- str_replace(paste(sample_name, truth_name, ifelse(use_roc_fdr, "fdr", ""), paste0(eventtype, collapse = "_"), sep="_"), "[ /]", "_")
all_ids <- c(truth_id, ids)
metadata <- LoadCachedMetadata(datadir)
metadata <- metadata %>% filter(Id %in% all_ids)
metadata$CX_REFERENCE_VCF <- list.files(datadir, pattern=paste0("^", truth_id, ".*.vcf$"))
callgr_na12878 <- .CachedLoadCallMatrixForIds(
		datadir=datadir,
		metadata=metadata,
		ids=all_ids,
		ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize,
		maxgap=maxgap, sizemargin=sizemargin, ignore.strand=ignore.strand,
		grtransform=lroptions$grtransform[[grtransformName]],
		grtransformName=grtransformName,
		nominalPosition=nominalPosition,
		eventtype=eventtype)

testdf_na12878 <- as.data.frame(callgr_na12878) %>%
	filter(!ignore.filtered) %>%
	filter(Id != truth_id) %>%
	mutate(tp=ifelse(fId00000000000000000000000000000003 >=0, "tp", "fp")) %>%
	filter_to_variables()

na12878_predictors <- testdf_na12878 %>% select(-tp, -svtype)
mprob_na12878 <- predict(chm13_rf, na12878_predictors, type="prob")

rocdf_na12878 <- as.data.frame(callgr_na12878) %>%
	filter(ignore.filtered) %>%
	filter(Id != truth_id) %>%
	mutate(tp=ifelse(fId00000000000000000000000000000003 >=0, "tp", "fp")) %>%
	filter_to_variables() %>%
	mutate(scoring="QUAL score (default filters)") %>%
	select(tp, scoring, QUAL) %>%
	rbind(testdf_na12878 %>% mutate(scoring="QUAL score (all calls)") %>%
		select(tp, scoring, QUAL)
	) %>%
	rbind(testdf_na12878 %>%
			select(tp) %>%
			mutate(scoring="Random forest", QUAL=mprob_na12878[,2])
	) %>%
	mutate(tp=tp == "tp") %>%
	mutate(fp=!tp) %>%
	group_by(scoring, QUAL) %>%
	summarise(fp=sum(fp), tp=sum(tp)) %>%
	group_by(scoring) %>%
	arrange(desc(QUAL)) %>%
	mutate(
		fp=cumsum(fp),
		tp=cumsum(tp)) %>%
	# each QUAL score is a point on the ROC plott
	group_by(scoring, QUAL) %>%
	summarise(tp=max(tp), fp=max(fp)) %>%
	# QUAL scores with the same number of tp calls can be merged on the ROC plot
	group_by(scoring, tp) %>%
	summarise(fp=max(fp), QUAL=min(QUAL)) %>%
	# subsample along tp and tp+fp axis
	group_by(scoring) %>%
	dplyr::slice(unique(c(
		1,
		findInterval(seq(0, max(tp), max(tp)/1000), tp),
		findInterval(seq(0, max(tp + fp), max(tp + fp)/1000), tp + fp),
		n()
	))) %>%
	mutate(is_endpoint = tp == max(tp)) %>%
	ungroup() %>%
	#TODO: fn and sens need eventCount using group_by(Id, CallSet, !!!groupingCols)
	mutate(
		precision=tp / (tp + fp),
		fdr=1-precision)
ggplot(rocdf_na12878) +
	aes(x=tp, y=precision, colour=scoring) +
	geom_line() +
	labs(title="NA12878 results based on random forest trained on CHM13\n(QUAL, svLen, insLen, HOMLEN, repeatClass, snp50bp)")


