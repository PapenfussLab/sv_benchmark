source("global.R")

plotdata <- NULL
# cache the results for all the knobs that can be turned in the UI
for (data in simoptions$ds) {
	for (maxgap in simoptions$maxgap) {
		for (ignore.strand in simoptions$ignore.strand) {
			for (sizemargin in simoptions$sizemargin) {
				for (ignore.duplicates in simoptions$ignore.duplicates) {
					for (ignore.interchromosomal in simoptions$ignore.interchromosomal) {
						for (mineventsize in simoptions$mineventsize) {
							#for (maxeventsize in simoptions$maxeventsize) {
								plotdata <- LoadPlotData(
									datadir=paste0("../data.", data),
									maxgap=maxgap,
									ignore.strand=ignore.strand,
									sizemargin=sizemargin,
									ignore.duplicates=ignore.duplicates,
									ignore.interchromosomal=ignore.interchromosomal,
									mineventsize=mineventsize,
									maxeventsize=NULL,
									vcftransform=simoptions$vcftransform,
									truthgr=NULL,
									existingCache=plotdata)
							#}
						}
					}
				}
			}
		}
	}
}
