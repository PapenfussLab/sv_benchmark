#!/usr/bin/env Rscript
library("optparse")
source("global.R")
source("libplot.R")

plotGraphs <- FALSE
loadAll <- TRUE
if (!interactive()) {
	# check command-line args and only process that data subset
	args <- commandArgs(trailingOnly=TRUE)
	simoptions$datadir <- simoptions$datadir[simoptions$datadir %in% args]
	lroptions$datadir <- lroptions$datadir[lroptions$datadir %in% args]

	parser <- OptionParser(option_list=list(
		make_option(c("-p", "--plot"), action="store_true", default=FALSE, help="Generate plots")
	))
	arguments = parse_args(parser, positional_arguments=TRUE)
	opt <- arguments$options
	args <- arguments$args
	plotGraphs <- opt$plot
	loadAll <- FALSE
}

plotdata <- NULL
# cache the results for all the knobs that can be turned in the UI
for (datadir in dataoptions$datadir[dataoptions$datadir %in% simoptions$datadir]) {
	for (maxgap in simoptions$maxgap) {
		for (ignore.strand in simoptions$ignore.strand) {
			for (sizemargin in simoptions$sizemargin) {
				for (ignore.duplicates in simoptions$ignore.duplicates) {
					for (ignore.interchromosomal in simoptions$ignore.interchromosomal) {
						for (requiredHits in simoptions$requiredHits) {
							for (mineventsize in simoptions$mineventsize) {
								transformName <- "PrimaryHumanOnly"
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
										grtransformName = transformName,
										grtransform = transform,
										requiredHits = requiredHits,
										truthbedpedir = NULL,
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
									# Event Size
									plotdf <- plotdata$dfs$callsByEventSize %>%
										filter(abs(svLen) >= mineventsize) %>% # hack to 48 shows up for 51bp min
										filter(StripCallerVersion(CX_CALLER, FALSE) %in% fulldatacallers) %>%
										inner_join(plotdata$dfs$mostSensitiveAligner) %>%
										mutate(
											CallSet=ifelse(CallSet=="High & Low confidence", "All Calls", CallSet),
											caller=StripCallerVersion(CX_CALLER, FALSE),
											eventtype=factor(eventtype, c("Deletion","Insertion","Inversion","Tandem Duplication","Breakpoint")))
									ggplot(plotdf) +
										aes_string(colour=paste0("as.factor(", simcolour, ")")) +
										#geom_label(data=expand.grid(caller=unique(plotdf$caller), eventtype=unique(plotdf$eventtype)), colour="grey", x=2**(16/2)-mineventsize, y=0.5, aes(label=caller)) +
										geom_line(size=0.5, aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, linetype=CallSet)) +
										scale_x_svlen +
										#coord_cartesian(xlim=c(min(abs(plotdf$svLen), max(abs(plotdf$svLen))))) +
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
											CallSet=ifelse(CallSet=="High & Low confidence", "All Calls", CallSet),
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
								}
							}
						}
					}
				}
			}
		}
	}
}
for (datadir in lroptions$datadir) {
	datapath = paste0(dataLocation, "data.", datadir)
	for (truthpath in lroptions$truthpath) {
		truthbedpedir = paste0(dataLocation, "input.", datadir, "/", truthpath)
		for (maxgap in lroptions$maxgap) {
			for (ignore.strand in lroptions$ignore.strand) {
				for (sizemargin in lroptions$sizemargin) {
					for (ignore.duplicates in lroptions$ignore.duplicates) {
						for (ignore.interchromosomal in lroptions$ignore.interchromosomal) {
							for (mineventsize in lroptions$mineventsize) {
								for (requiredHits in lroptions$requiredHits) {
									for (transformName in names(lroptions$grtransform)) {
										transform <- lroptions$grtransform[[transformName]]
									  for (et in eventtypes) {
											plotdata <- LoadPlotData(
													datadir = datapath,
													maxgap = maxgap,
													ignore.strand = ignore.strand,
													sizemargin = sizemargin,
													ignore.duplicates = ignore.duplicates,
													ignore.interchromosomal = ignore.interchromosomal,
													mineventsize = mineventsize,
													maxeventsize = lroptions$maxeventsize,
													grtransformName = transformName,
													grtransform = transform,
													requiredHits = requiredHits,
													truthbedpedir = truthbedpedir,
													eventtypes=et,
													existingCache = plotdata,
													loadFromCacheOnly=FALSE,
													loadAll=loadAll)
											if (plotGraphs) {
												filenamePrefix <- paste(datadir, et, maxgap, ignore.strand, sizemargin, ignore.duplicates, ignore.interchromosomal, requiredHits, mineventsize, transformName, str_replace(truthpath, "/", "_"), sep="_")
												####
												# Precision-recall
												plotdf <- plotdata$dfs$roc %>%
													filter(StripCallerVersion(CX_CALLER, FALSE) %in% fulldatacallers) %>%
													inner_join(plotdata$dfs$mostSensitiveAligner) %>%
													mutate(
														CallSet=ifelse(CallSet=="High & Low confidence", "All Calls", CallSet),
														caller=StripCallerVersion(CX_CALLER, FALSE)) %>%
													rbind(data.frame(Id=paste0(fulldatacallers, "_placeholder"), CallSet="High confidence only", tp=0, events=1, fp=0, fn=0, QUAL=-1, precision=0, fdr=0, sens=0,
																					 CX_ALIGNER=NA,CX_ALIGNER_MODE="",CX_MULTIMAPPING_LOCATIONS=NA,CX_CALLER=fulldatacallers,CX_READ_LENGTH=NA,CX_READ_DEPTH=NA,CX_READ_FRAGMENT_LENGTH=NA,CX_REFERENCE_VCF=NA,
																					 caller=fulldatacallers))
												ggplot(plotdf %>% arrange(desc(QUAL))) +
													aes(group = paste(Id, CallSet), y = precision, x = tp, colour=caller, linetype=CallSet) +
													geom_line(size=1) +
													scale_colour_brewer(palette = "Paired") +
													labs(title = paste(datadir, names(eventtypes[eventtypes==et])), y = "Precision", x = "Recall (true positive count)")
												saveplot(paste0("precrecall_", filenamePrefix), scale=2, height=12, width=12, units="cm")
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

