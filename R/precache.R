#!/usr/bin/env Rscript
source("global.R")

if (!interactive()) {
	# check command-line args and only process that data subset
	args <- commandArgs(trailingOnly=TRUE)
	simoptions$datadir <- simoptions$datadir[simoptions$datadir %in% args]
	lroptions$datadir <- lroptions$datadir[lroptions$datadir %in% args]
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
								transform <- simoptions$grtransform
								plotdata <- LoadPlotData(
										datadir = paste0(dataLocation, "data.", datadir),
										maxgap = maxgap,
										ignore.strand = ignore.strand,
										sizemargin = sizemargin,
										ignore.duplicates = ignore.duplicates,
										ignore.interchromosomal = ignore.interchromosomal,
										mineventsize = mineventsize,
										maxeventsize = simoptions$maxeventsize,
										grtransform = transform,
										requiredHits = requiredHits,
										truthgr = NULL,
										eventtypes=NULL,
										existingCache = plotdata)
							}
						}
					}
				}
			}
		}
	}
}
for (datadir in lroptions$datadir) {
	truthgr <- LoadLongReadTruthgr(paste0(dataLocation, "input.", datadir, "/longread"))
	for (maxgap in lroptions$maxgap) {
		for (ignore.strand in lroptions$ignore.strand) {
			for (sizemargin in lroptions$sizemargin) {
				for (ignore.duplicates in lroptions$ignore.duplicates) {
					for (ignore.interchromosomal in lroptions$ignore.interchromosomal) {
						for (mineventsize in lroptions$mineventsize) {
							for (requiredHits in lroptions$requiredHits) {
								for (transform in lroptions$grtransform) {
								  # UI now drop-down so only need for each event type
								  # all possible event type combinations
									#for (i in 1:(2**length(eventtypes)-1)) { 
								  #et <- eventtypes[as.logical(intToBits(i))[1:length(eventtypes)]]
								  for (et in eventtypes) {
										plotdata <- LoadPlotData(
												datadir = paste0(dataLocation, "data.", datadir),
												maxgap = maxgap,
												ignore.strand = ignore.strand,
												sizemargin = sizemargin,
												ignore.duplicates = ignore.duplicates,
												ignore.interchromosomal = ignore.interchromosomal,
												mineventsize = mineventsize,
												maxeventsize = lroptions$maxeventsize,
												grtransform = transform,
												requiredHits = requiredHits,
												truthgr = truthgr,
												eventtypes=et,
												existingCache = plotdata)
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

