source("global.R")
library(shiny)
library(ggplot2)
library(dplyr)

data <- NULL

shinyServer(function(input, output) {
  output$mainPlot <- renderPlot({
    mineventsize <- NULL
    if (!input$smallevents) mineventsize <- 51
  	data <- LoadPlotData(
  		datadir=paste0("../data.", input$data),
  		maxgap=200,
  		ignore.strand=TRUE,
  		sizemargin=0.25,
  		ignore.duplicates=TRUE,
  		ignore.interchromosomal=TRUE,
  		mineventsize=mineventsize,
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
  	mostSens <- paste(data$dfs$mostSensitiveAligner$Id, data$dfs$mostSensitiveAligner$CallSet)
		es <- data$dfs$callsByEventSize
		es <- es %>%
			filter(
				CX_REFERENCE_VCF_VARIANTS %in% input$eventtype &
				CX_READ_LENGTH %in% as.numeric(input$readlength) &
				CX_READ_DEPTH %in% as.numeric(input$depth) &
				CX_READ_FRAGMENT_LENGTH %in% as.numeric(input$fragsize))
  	# aligner filters
  	alignerIdCallSet <- es %>%
  			select(Id, CallSet, CX_ALIGNER) %>%
  			filter(CX_ALIGNER %in% input$aligner | (is.na(CX_ALIGNER) & "" %in% input$aligner)) %>%
  			select(Id, CallSet) %>%
  			rbind(data$dfs$mostSensitiveAligner[rep("best" %in% input$aligner, nrow(data$dfs$mostSensitiveAligner)),]) %>%
  			distinct()
		es <- es %>% inner_join(alignerIdCallSet)
  	es <- es %>% mutate(
  	  caller=StripCallerVersion(CX_CALLER),
  	  aligner=CX_ALIGNER %na% "N/A"           )
		p <- ggplot(es) +
			aes(group=paste(Id, CallSet), x=abs(svLen), y=sens) +
  		aes_string(linetype=paste0("as.factor(", input$linetype, ")"), colour=input$colour) +
			geom_line(size=0.5) +
			scale_x_svlen +
			facet_grid(caller ~ CX_READ_DEPTH + CX_READ_LENGTH + CX_READ_FRAGMENT_LENGTH) +
			labs(title="", y="Sensitivity", x="Event size",
				linetype=names(facets)[facets==input$linetype],
				colour=names(facets)[facets==input$colour])
  	return(p)
  })
})
