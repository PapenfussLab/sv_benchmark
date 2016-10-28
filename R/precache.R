source("global.R")

plotdata <- NULL
# cache the results for all the knobs that can be turned in the UI
for (datadir in simoptions$datadir) {
	for (maxgap in simoptions$maxgap) {
		for (ignore.strand in simoptions$ignore.strand) {
			for (sizemargin in simoptions$sizemargin) {
				for (ignore.duplicates in simoptions$ignore.duplicates) {
					for (ignore.interchromosomal in simoptions$ignore.interchromosomal) {
						for (requiredHits in lroptions$requiredHits) {
							for (mineventsize in simoptions$mineventsize) {
								plotdata <- LoadPlotData(
										datadir = paste0(dataLocation, "data.", datadir),
										maxgap = maxgap,
										ignore.strand = ignore.strand,
										sizemargin = sizemargin,
										ignore.duplicates = ignore.duplicates,
										ignore.interchromosomal = ignore.interchromosomal,
										mineventsize = mineventsize,
										maxeventsize = simoptions$maxeventsize,
										vcftransform = simoptions$vcftransform,
										requiredHits = requiredHits,
										truthgr = NULL,
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
	truthgr <- LoadLongReadTruthgr(paste0(dataLocation, "input.", datadir))
	for (maxgap in lroptions$maxgap) {
		for (ignore.strand in lroptions$ignore.strand) {
			for (sizemargin in lroptions$sizemargin) {
				for (ignore.duplicates in lroptions$ignore.duplicates) {
					for (ignore.interchromosomal in lroptions$ignore.interchromosomal) {
						for (mineventsize in lroptions$mineventsize) {
							for (requiredHits in lroptions$requiredHits) {
								plotdata <- LoadPlotData(
										datadir = paste0(dataLocation, "data.", datadir),
										maxgap = maxgap,
										ignore.strand = ignore.strand,
										sizemargin = sizemargin,
										ignore.duplicates = ignore.duplicates,
										ignore.interchromosomal = ignore.interchromosomal,
										mineventsize = mineventsize,
										maxeventsize = lroptions$maxeventsize,
										vcftransform = lroptions$vcftransform,
										requiredHits = requiredHits,
										truthgr = NULL,
										existingCache = plotdata)
							}
						}
					}
				}
			}
		}
	}
}
