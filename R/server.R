source("sv_benchmark.R")
library(shiny)
library(ggplot2)
library(dplyr)
mvcfs <- NULL
callset <- NULL
gdf <- NULL
shinyServer(function(input, output) {
  output$mainPlot <- renderPlot({
  	mvcfs <- LoadVCFs("../data.aligner", result=mvcfs)
		callset <- LoadSimulatedCallSets(mvcfs, maxgap=200, ignore.strand=TRUE, sizemargin=0.25, result=callset)
  	gdf <- LoadGraphDataFrames(callset, ignore.duplicates=TRUE, ignore.interchromosomal=TRUE, mineventsize=NULL, maxeventsize=NULL, result=gdf)
  	ggplot(gdf$callsByEventSize) +
  		aes(group=paste(Id, CallSet), x=abs(svLen), y=sens, color=eventtype, linetype=CallSet) +
		  geom_line(size=0.5) +
		  scale_x_svlen +
		  facet_grid(CX_CALLER ~ CX_READ_DEPTH, switch="y") +
		  labs(title="", y="Sensitivity", x="Event size", color="Event Type", linetype="Call set")
  })
})
