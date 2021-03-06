#!/usr/bin/env Rscript
#library(parallel)
library("optparse")
source("global.R")
source("libplot.R")
source("server.R")

# check command-line args and only process that data subset
args <- commandArgs(trailingOnly=TRUE)
parser <- OptionParser(option_list=list(
	make_option(c("--qsub"), action="store_true", default=FALSE, help="Output qsub submission scripts."),
	make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate plots"),
	make_option(c("--datadir"), type="character"),
	make_option(c("--maxgap"), type="integer"),
	make_option(c("--ignore.strand"), type="logical"),
	make_option(c("--sizemargin"), type="integer"),
	make_option(c("--ignore.duplicates"), type="logical"),
	make_option(c("--ignore.interchromosomal"), type="logical"),
	make_option(c("--mineventsize"), type="integer"),
	make_option(c("--maxeventsize"), type="integer"),
	make_option(c("--grtransformName"), type="character"),
	make_option(c("--requiredHits"), type="integer"),
	make_option(c("--truthbedpedir"), type="character"),
	make_option(c("--mintruthbedpescore"), type="integer")
))
arguments = parse_args(parser)
plotGraphs <- arguments$plot
loadAll <- FALSE

subsetToArgs <- function(df) {
	for (argName in names(arguments)) {
		if (!is.null(df[[argName]])) {
			df <- df[df[[argName]] %in% arguments[[argName]],]
		}
	}
	return(df)
}
# To generate the supp figures in the paper:
simoptions$maxgap <- 200
simoptions$sizemargin <- 0.25
simoptions$ignore.strand <- TRUE
simoptions$mineventsize <- 48
simoptions$ignore.interchromosomal
simoptions$ignore.duplicates <- TRUE
rootdir <- "W:/i/"
dataLocation <- "W:/i/"
dataoptions$datadir <- simoptions$datadir


plotdata <- NULL
simcache <- function(datadir, maxgap, ignore.strand, sizemargin, ignore.duplicates, ignore.interchromosomal, requiredHits, mineventsize) {
	write(paste0("Processing ",
		",datadir=", datadir,
		",maxgap=", maxgap,
		",ignore.strand=", ignore.strand,
		"sizemargin=", sizemargin,
		"ignore.duplicates=", ignore.duplicates,
		"ignore.interchromosomal=", ignore.interchromosomal,
		"requiredHits=", requiredHits,
		"mineventsize=", mineventsize))
	grtransformName <- "PrimaryHumanOnly"
	transform <- simoptions$grtransform[["PrimaryHumanOnly"]]
	plotdata <- LoadPlotData(
			datadir = paste0(dataLocation, "data.", datadir),
			maxgap = maxgap,
			ignore.strand = ignore.strand,
			sizemargin = sizemargin,
			ignore.duplicates = ignore.duplicates,
			ignore.interchromosomal = ignore.interchromosomal,
			mineventsize = mineventsize,
			maxeventsize = simoptions$maxeventsize,
			grtransformName = grtransformName,
			grtransform = transform,
			requiredHits = requiredHits,
			truthbedpedir = NULL,
			mintruthbedpescore=NULL,
			eventtypes=NULL,
			existingCache = plotdata,
			loadFromCacheOnly=FALSE,
			loadAll=loadAll)
	if (plotGraphs) {
		filenamePrefix <- paste(datadir, maxgap, ignore.strand, sizemargin, ignore.duplicates, ignore.interchromosomal, requiredHits, mineventsize, sep="_")
		simcolour <-list(
			"rl"="CX_READ_LENGTH",
			"rd"="CX_READ_DEPTH",
			"fs"="CX_READ_FRAGMENT_LENGTH")[[datadir]]
		simlinetype <- "CallSet"
		####
		# Breakend error distribution
		plotdata$dfs$bpErrorDistribution %>%
		  filter(CallSet==ALL_CALLS) %>%
		  mutate(fill=if_else(nominalPosition, "Nominal position", "Incorporating caller\nconfidence interval and\nreported microhomology")) %>%
		  mutate(
		    caller=StripCallerVersion(CX_CALLER, FALSE),
		    eventtype=factor(eventtype, c("Deletion","Insertion","Inversion","Tandem Duplication","Breakpoint"))) %>%
		ggplot() +
		  aes(x=bperror, y=rate, fill=fill) +
		  #geom_bar(stat="identity") +
		  geom_bar(stat="identity", data=plotdf %>% filter(nominalPosition==FALSE), alpha=0.5) +
		  geom_bar(stat="identity", data=plotdf %>% filter(nominalPosition==TRUE), alpha=0.5) +
		  scale_fill_brewer(type="qual", palette="Dark2", name="Error margin") +
		  facet_wrap(eventtype ~ caller) +
		  labs(title="Error in called position", x="Error (base pairs)", y="Portion of true positive calls")
		saveplot(paste0("sim_error_", filenamePrefix), scale=2, height=18, width=12, units="cm")
		####
		# Event Size
		plotdf <- plotdata$dfs$callsByEventSize %>%
			filter(abs(svLen) >= mineventsize) %>% # hack to 48 shows up for 51bp min
			filter(StripCallerVersion(CX_CALLER, FALSE) %in% fulldatacallers) %>%
			inner_join(plotdata$dfs$mostSensitiveAligner) %>%
			mutate(
				CallSet=ifelse(CallSet==ALL_CALLS, "All Calls", CallSet),
				caller=StripCallerVersion(CX_CALLER, FALSE),
				eventtype=factor(eventtype, c("Deletion","Insertion","Inversion","Tandem Duplication","Breakpoint")))
		ggplot(plotdf) +
			aes_string(colour=paste0("as.factor(", simcolour, ")")) +
			#geom_label(data=expand.grid(caller=unique(plotdf$caller), eventtype=unique(plotdf$eventtype)), colour="grey", x=2**(16/2)-mineventsize, y=0.5, aes(label=caller)) +
			geom_line(size=0.5, aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, linetype=CallSet)) +
			scale_x_svlen +
			coord_cartesian(xlim=c(50, 2**16)) + # start at 50bp as per reviewer request
			#scale_colour_brewer(palette="Dark2") +
			# simfacets for the fields not displayed in linetype or colour
			facet_grid(caller ~ eventtype) +
			labs(title="", y="Sensitivity", x="Event size", linetype="Call Set", colour=names(simfacets)[simfacets==simcolour]) +
			theme(axis.text.x = element_text(angle=-90, vjust=0.5, hjust=0))
		saveplot(paste0("simeventsize_", filenamePrefix), scale=2, height=18, width=12, units="cm")
		####
		# ROC
		plotdf <- plotdata$dfs$roc %>%
			filter(StripCallerVersion(CX_CALLER, FALSE) %in% fulldatacallers) %>%
			filter(!(CX_REFERENCE_VCF_VARIANTS %in% "hetBP_SINE")) %>%
			inner_join(plotdata$dfs$mostSensitiveAligner) %>%
			mutate(
				CallSet=ifelse(CallSet==ALL_CALLS, "All Calls", CallSet),
				caller=StripCallerVersion(CX_CALLER, FALSE),
				eventtype=factor(eventtype, c("Deletion","Insertion","Inversion","Tandem Duplication","Breakpoint")))
		ggplot(plotdf) +
			aes(group=paste(Id, CallSet), y=sens, x=fp+1, linetype=CallSet) +
			aes_string(colour=paste0("as.factor(", simcolour, ")")) +
			#geom_label(data=expand.grid(caller=unique(plotdf$caller), eventtype=unique(plotdf$eventtype)), colour="grey", x=2**(16/2)-mineventsize, y=0.5, aes(label=caller)) +
			geom_line() +
			scale_x_log_fp +
			coord_cartesian(xlim=c(1, max(plotdf$fp))) +
			#coord_cartesian(xlim=c(min(abs(plotdf$svLen), max(abs(plotdf$svLen))))) +
			#scale_colour_brewer(palette="Dark2") +
			# simfacets for the fields not displayed in linetype or colour
			facet_grid(caller ~ eventtype) +
			labs(title="", y="Sensitivity", x="False Positives", linetype="Call Set", colour=names(simfacets)[simfacets==simcolour])
		saveplot(paste0("simroc_", filenamePrefix), scale=2, height=18, width=12, units="cm")
		browser()
	}
}
simparam = expand.grid(
	mineventsize=simoptions$mineventsize,
	ignore.interchromosomal=simoptions$ignore.interchromosomal,
	ignore.duplicates=simoptions$ignore.duplicates,
	maxgap=simoptions$maxgap,
	ignore.strand=simoptions$ignore.strand,
	sizemargin=simoptions$sizemargin,
	requiredHits=simoptions$requiredHits,
	datadir=dataoptions$datadir[dataoptions$datadir %in% simoptions$datadir],
	stringsAsFactors=FALSE)
simparam <-simparam[sample(nrow(simparam)),] # randomise row ordering
simparam <- subsetToArgs(simparam)
simcachedf <- function(pass, df) {
	if (nrow(df) > 0) {
		if (arguments$qsub) {
			write("# cache blocks that can be executed in parallel:", stdout())
			for (i in 1:nrow(df)) {
				write(paste(
					"echo Rscript precache.R",
					ifelse(arguments$plot, "--plot", ""),
					"--datadir", df$datadir[i],
					"--maxgap", df$maxgap[i],
					"--ignore.strand", df$ignore.strand[i],
					"--sizemargin", df$sizemargin[i],
					"--ignore.duplicates", df$ignore.duplicates[i],
					"--ignore.interchromosomal", df$ignore.interchromosomal[i],
					"--mineventsize", df$mineventsize[i],
					"--requiredHits", df$requiredHits[i],
					"|",
					"qsub",
					"-V",
					"-D", "$PWD",
					"-N", paste0(pass, "-", df$datadir[i])
					), stdout())
			}
		} else {
			lapply(seq_len(nrow(df)), function(i) {
				simcache(df$datadir[i], df$maxgap[i], df$ignore.strand[i], df$sizemargin[i], df$ignore.duplicates[i], df$ignore.interchromosomal[i], df$requiredHits[i], df$mineventsize[i])
				return(0)
			})#, mc.preschedule=FALSE)
		}
	}
}
# Since R.cache is not thread safe, we need to ensure we don't run two operations
# which both write to the same cached file at the same time
simcachedf("p1", simparam %>% distinct(datadir, .keep_all=TRUE))
simcachedf("p2", simparam %>% distinct(datadir, requiredHits, sizemargin, ignore.strand, maxgap, .keep_all=TRUE))
simcachedf("p3", simparam)

lrcache <- function(datadir, truthbedpedir, mintruthbedpescore, maxgap, ignore.strand, sizemargin, requiredHits, grtransformName, ignore.duplicates, ignore.interchromosomal, mineventsize) {
	input <- list(
		datasettype="lr",
		datadir=datadir,
		lrmaxgap=maxgap,
		lrignore.strand=ignore.strand,
		lrsizemargin=sizemargin,
		lrignore.duplicates=ignore.duplicates,
		ignore.interchromosomal=ignore.interchromosomal,
		mineventsize=mineventsize,
		maxeventsize=NULL,
		lrgrtransformName=grtransformName,
		truthbedpedir=truthbedpedir,
		mintruthbedpescore=mintruthbedpescore,
		lrTruthSet=ifelse(is.null(truthbedpedir), "Published Calls", "Long Reads"),
		lraligner="best",
		lrcaller=fulldatacallers
		)
	datapath = paste0(dataLocation, "data.", datadir)
	truthbedpedir = paste0(dataLocation, "input.", datadir, "/", truthbedpedir)
	transform <- lroptions$grtransform[[grtransformName]]
  for (et in eventtypes) {
  	input$eventtypes <- et
		plotdata <- LoadPlotData(
				datadir = datapath,
				maxgap = maxgap,
				ignore.strand = ignore.strand,
				sizemargin = sizemargin,
				ignore.duplicates = ignore.duplicates,
				ignore.interchromosomal = ignore.interchromosomal,
				mineventsize = mineventsize,
				maxeventsize = lroptions$maxeventsize,
				grtransformName = grtransformName,
				grtransform = transform,
				requiredHits = requiredHits,
				truthbedpedir = truthbedpedir,
				mintruthbedpescore=mintruthbedpescore,
				eventtypes=et,
				existingCache=plotdata,
				loadFromCacheOnly=FALSE,
				loadAll=loadAll)
		if (plotGraphs) {
			filenamePrefix <- str_replace_all(paste(datadir, et, maxgap, ignore.strand, sizemargin, ignore.duplicates, ignore.interchromosomal, requiredHits, mineventsize, grtransformName, truthbedpedir, sep="_"), "/", "_")
			saveplot(paste0(filenamePrefix, "_bpErrorDistribution"), scale=2, height=12, width=12, units="cm", plot=doPlot(input, "bpErrorDistribution", plotCallPositionErrorDistribution, debugLabel="bpErrorDistributionPlot", data=plotdata))
			saveplot(paste0(filenamePrefix, "_precrecall"), scale=2, height=12, width=12, units="cm", plot=doPlot(input, "roc", plotPrecRecall, debugLabel="lrPrecRecallPlot", data=plotdata))
			####
			# Precision-recall
			#plotdf <- plotdata$dfs$roc %>%
			#	filter(StripCallerVersion(CX_CALLER, FALSE) %in% fulldatacallers) %>%
			#	inner_join(plotdata$dfs$mostSensitiveAligner) %>%
			#	mutate(
			#		CallSet=ifelse(CallSet==ALL_CALLS, "All Calls", CallSet),
			#		caller=StripCallerVersion(CX_CALLER, FALSE)) %>%
			#	rbind(data.frame(Id=paste0(fulldatacallers, "_placeholder"), CallSet=PASS_CALLS, tp=0, events=1, fp=0, fn=0, QUAL=-1, precision=0, fdr=0, sens=0,
			#									 CX_ALIGNER=NA,CX_ALIGNER_MODE="",CX_MULTIMAPPING_LOCATIONS=NA,CX_CALLER=fulldatacallers,CX_READ_LENGTH=NA,CX_READ_DEPTH=NA,CX_READ_FRAGMENT_LENGTH=NA,CX_REFERENCE_VCF=NA,
			#									 caller=fulldatacallers))
			#ggplot(plotdf %>% arrange(desc(QUAL))) +
			#	aes(group = paste(Id, CallSet), y = precision, x = tp, colour=caller, linetype=CallSet) +
			#	geom_line(size=1) +
			#	scale_colour_brewer(palette = "Paired") +
			#	labs(title = paste(datadir, names(eventtypes[eventtypes==et])), y = "Precision", x = "Recall (true positive count)")
			#saveplot(paste0("precrecall_", filenamePrefix), scale=2, height=12, width=12, units="cm")
		}
	}
}

lrparam = expand.grid(
	# Rcache/LoadPlotData
	ignore.duplicates=lroptions$ignore.duplicates,
	ignore.interchromosomal=lroptions$ignore.interchromosomal,
	mineventsize=lroptions$mineventsize,
	# .Rcache/ScoreVariantsFromTruthVCF
	maxgap=lroptions$maxgap,
	ignore.strand=lroptions$ignore.strand,
	sizemargin=lroptions$sizemargin,
	requiredHits=lroptions$requiredHits,
	grtransformName=names(lroptions$grtransform),
	# .Rcache/import.sv.bedpe
	truthbedpedir=lroptions$truthpath,
	mintruthbedpescore=lroptions$mintruthscore,
	# .Rcache/LoadMinimalSVs
	datadir=lroptions$datadir,
	stringsAsFactors=FALSE
	)
lrparam <-lrparam[sample(nrow(lrparam)),] # randomise row ordering
lrparam <- subsetToArgs(lrparam)
lrcachedf <- function(pass, df) {
	if (nrow(df) > 0) {
		if (arguments$qsub) {
			write("# cache blocks that can be executed in parallel:", stdout())
			for (i in 1:nrow(df)) {
				write(paste(
					"echo Rscript precache.R",
					ifelse(arguments$plot, "--plot", ""),
					"--datadir", df$datadir[i],
					"--maxgap", df$maxgap[i],
					"--ignore.strand", df$ignore.strand[i],
					"--sizemargin", df$sizemargin[i],
					"--ignore.duplicates", df$ignore.duplicates[i],
					"--ignore.interchromosomal", df$ignore.interchromosomal[i],
					"--mineventsize", df$mineventsize[i],
					"--grtransformName", df$grtransformName[i],
					"--requiredHits", df$requiredHits[i],
					"--truthbedpedir", df$truthbedpedir[i],
					"--mintruthbedpescore", df$mintruthbedpescore[i],
					"|",
					"qsub",
					"-V",
					"-D", "$PWD",
					"-N", paste0(pass, "-", df$datadir[i])
				),stdout())
			}
		} else {
			lapply(seq_len(nrow(df)), function(i) {
				# long read
				lrcache(df$datadir[i], df$truthbedpedir[i], df$mintruthbedpescore[i], df$maxgap[i], df$ignore.strand[i], df$sizemargin[i], df$requiredHits[i], df$grtransformName[i], df$ignore.duplicates[i], df$ignore.interchromosomal[i], df$mineventsize[i])
				# default call set
				lrcache(df$datadir[i], NULL, NULL, df$maxgap[i], df$ignore.strand[i], df$sizemargin[i], 1, df$grtransformName[i], df$ignore.duplicates[i], df$ignore.interchromosomal[i], df$mineventsize[i])
				return(0)
			})#, mc.preschedule=FALSE)
		}
	}
}
lrcachedf("p1", lrparam %>% distinct(datadir, .keep_all=TRUE))
lrcachedf("p2", lrparam %>% distinct(datadir, truthbedpedir, mintruthbedpescore, .keep_all=TRUE))
lrcachedf("p3", lrparam %>% distinct(datadir, truthbedpedir, mintruthbedpescore, maxgap, ignore.strand, sizemargin, requiredHits, grtransformName, .keep_all=TRUE))
lrcachedf("p4", lrparam)
