source("global.R")
library(shiny)
library(ggplot2)
library(dplyr)

data <- NULL
RefreshSimData <- function(input, olddata) {
	if (!input$smallevents) {
		mineventsize <- 51
	} else {
		mineventsize <- 0
	}
	LoadPlotData(
  		datadir=paste0("../data.", input$data),
  		maxgap=simoptions$maxgap,
  		ignore.strand=simoptions$ignore.strand,
  		sizemargin=simoptions$sizemargin,
  		ignore.duplicates=simoptions$ignore.duplicates,
  		ignore.interchromosomal=simoptions$ignore.interchromosomal,
  		mineventsize=mineventsize,
  		maxeventsize=simoptions$maxeventsize,
  		vcftransform=simoptions$vcftransform,
			truthgr=NULL,
			existingCache=olddata)
}
PrettyFormatPlotdf <- function(input, data, plotdf) {
	plotdf <- plotdf %>%
			filter(
				StripCallerVersion(CX_CALLER, FALSE) %in% input$caller &
				CallSet %in% input$callset &
				CX_REFERENCE_VCF_VARIANTS %in% input$eventtype &
				CX_READ_LENGTH %in% as.numeric(input$readlength) &
				CX_READ_DEPTH %in% as.numeric(input$depth) &
				CX_READ_FRAGMENT_LENGTH %in% as.numeric(input$fragsize))
  	# aligner filters
  	alignerIdCallSet <- plotdf %>%
  			select(Id, CallSet, CX_ALIGNER) %>%
  			filter(CX_ALIGNER %in% input$aligner | (is.na(CX_ALIGNER) & "" %in% input$aligner)) %>%
  			select(Id, CallSet) %>%
  			rbind(data$dfs$mostSensitiveAligner[rep("best" %in% input$aligner, nrow(data$dfs$mostSensitiveAligner)),]) %>%
  			distinct()
		plotdf <- plotdf %>% inner_join(alignerIdCallSet)
  	plotdf <- plotdf %>% mutate(
  	  caller=StripCallerVersion(CX_CALLER, FALSE),
  	  aligner=CX_ALIGNER %na% rep("N/A", nrow(plotdf)))
	return(plotdf)
}

shinyServer(function(input, output) {
	output$simulationControls <- renderUI({
		return(span(
			checkboxGroupInput("callset", "Call Set",
					c("High confidence only", "High & Low confidence"),
					"High & Low confidence"),
			selectInput("eventtype", "Event Type",
					withnames(sort(unique(md[[input$data]]$CX_REFERENCE_VCF_VARIANTS)), PrettyVariants(sort(unique(md[[input$data]]$CX_REFERENCE_VCF_VARIANTS)))),
					"hetDEL"),
				checkboxGroupInput("readlength", "Read Length (paired-end)",
					sort(unique(md[[input$data]]$CX_READ_LENGTH)),
					sort(unique(md[[input$data]]$CX_READ_LENGTH))),
				checkboxGroupInput("depth", "Read Depth (coverage)",
					sort(unique(md[[input$data]]$CX_READ_DEPTH)),
					sort(unique(md[[input$data]]$CX_READ_DEPTH))),
				checkboxGroupInput("fragsize", "Mean Library Fragment Size",
					sort(unique(md[[input$data]]$CX_READ_FRAGMENT_LENGTH)),
					sort(unique(md[[input$data]]$CX_READ_FRAGMENT_LENGTH))),
				checkboxGroupInput("aligner", "Aligner",
					PrettyAligner(input$data),
					"best"),
				checkboxGroupInput("caller", "Software",
					sort(as.character(unique(StripCallerVersion(md[[input$data]]$CX_CALLER)))),
					sort(as.character(unique(StripCallerVersion(md[[input$data]]$CX_CALLER)))))
			))
	})
  output$eventSizePlot <- renderPlot({
		data <- RefreshSimData(input, data)
  	plotdf <- PrettyFormatPlotdf(input, data, data$dfs$callsByEventSize)
  	if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf) +
			aes(group=paste(Id, CallSet), x=abs(svLen), y=sens) +
  		aes_string(linetype=paste0("as.factor(", input$linetype, ")"), colour=paste0("as.factor(", input$colour, ")")) +
			geom_line(size=0.5) +
			scale_x_svlen +
  		# facets for the fields not displayed in linetype or colour
			facet_grid(eval(parse(text=paste("caller ~ ", paste(facets[!(facets %in% c(input$linetype, input$colour))], collapse=" + "))))) +
			labs(title="", y="Sensitivity", x="Event size",
				linetype=names(facets)[facets==input$linetype],
				colour=names(facets)[facets==input$colour])
  	return(p)
  })
	output$rocPlot <- renderPlot({
		data <- RefreshSimData(input, data)
  	plotdf <- PrettyFormatPlotdf(input, data, data$dfs$roc)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
			aes(group=paste(Id, CallSet), y=sens, x=fp+1) +
  		aes_string(linetype=paste0("as.factor(", input$linetype, ")"), colour=paste0("as.factor(", input$colour, ")")) +
			geom_line() +
			scale_x_log_fp +
			facet_grid(eval(parse(text=paste("caller ~ ", paste(facets[!(facets %in% c(input$linetype, input$colour))], collapse=" + "))))) +
			labs(title="", y="Sensitivity", x="False Positives",
				linetype=names(facets)[facets==input$linetype],
				colour=names(facets)[facets==input$colour])
  	return(p)
  })
})
