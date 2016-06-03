source("sv_benchmark.R")
source("shinyCache.R")
source("libplot.R")
library(shiny)
library(ggplot2)
library(dplyr)

data <- NULL

shinyServer(function(input, output) {
  output$mainPlot <- renderPlot({
  	data <- LoadPlotData(
  		datadir="../data.test",
  		maxgap=200,
  		ignore.strand=TRUE,
  		sizemargin=0.25,
  		ignore.duplicates=TRUE,
  		ignore.interchromosomal=TRUE,
  		mineventsize=NULL,
  		maxeventsize=NULL,
  		vcftransform=function(id, metadata, vcfs) {
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
			},
  		truthgr=NULL,
  		existingCache=data)
  	p <- ggplot(data$dfs$callsByEventSize) +
  		aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, color=eventtype, linetype=CallSet) +
		  geom_line(size=0.5) +
		  scale_x_svlen +
		  facet_grid(CX_CALLER ~ CX_READ_DEPTH, switch="y") +
		  labs(title="", y="Sensitivity", x="Event size", color="Event Type", linetype="Call set")
  	return(p)
  })
})
