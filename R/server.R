source("global.R")
library(shiny)
library(ggplot2)
library(dplyr)

cachedsimdata <- NULL
cachedlrdata <- NULL
RefreshSimData <- function(input, olddata) {
	if (!input$smallevents) {
		mineventsize <- 51
	} else {
		mineventsize <- 0
	}
	LoadPlotData(
		datadir=paste0(dataLocation, "data.", input$simdatadir),
		maxgap=simoptions$maxgap,
		ignore.strand=simoptions$ignore.strand,
		sizemargin=simoptions$sizemargin,
		ignore.duplicates=simoptions$ignore.duplicates,
		ignore.interchromosomal=simoptions$ignore.interchromosomal,
		mineventsize=mineventsize,
		maxeventsize = simoptions$maxeventsize,
		requiredHits=1,
		vcftransform=simoptions$vcftransform,
		truthgr=NULL,
		existingCache = olddata)
}
RefreshlrData <- function(input, olddata) {
    LoadPlotData(
        datadir = paste0(dataLocation, "data.", input$lrdatadir),
        maxgap = lroptions$maxgap,
        ignore.strand = lroptions$ignore.strand,
        sizemargin = lroptions$sizemargin,
        ignore.duplicates = lroptions$ignore.duplicates,
        ignore.interchromosomal = lroptions$ignore.interchromosomal,
        mineventsize = lroptions$mineventsize,
        maxeventsize = lroptions$maxeventsize,
        requiredHits = lroptions$requiredHits,
				vcftransform = lroptions$vcftransform,
				truthgr = LoadLongReadTruthgr(paste0(dataLocation, "input.", input$lrdatadir, "/longread")),
        existingCache = olddata)
}
PrettyFormatSimPlotdf <- function(input, simdata, plotdf) {
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
				rbind(simdata$dfs$mostSensitiveAligner[rep("best" %in% input$aligner, nrow(simdata$dfs$mostSensitiveAligner)),]) %>%
				distinct()
		plotdf <- plotdf %>% inner_join(alignerIdCallSet)
		plotdf <- plotdf %>% mutate(
			caller=StripCallerVersion(CX_CALLER, FALSE),
			aligner=CX_ALIGNER %na% rep("N/A", nrow(plotdf)))
	return(plotdf)
}

# server function (must be last in file)
function(input, output, session) {
	setBookmarkExclude(c("bookmark"))
	observe({
		# Trigger this observer every time an input changes
		reactiveValuesToList(input)
		session$doBookmark()
	})
	onBookmarked(function(url) {
		updateQueryString(url)
	})
	output$simControls <- renderUI({
		return(span(
			selectInput("eventtype", "Event Type",
					withnames(sort(unique(md[[input$simdatadir]]$CX_REFERENCE_VCF_VARIANTS)), PrettyVariants(sort(unique(md[[input$simdatadir]]$CX_REFERENCE_VCF_VARIANTS)))),
					"hetDEL"),
				checkboxGroupInput("readlength", "Read Length (paired-end)",
					sort(unique(md[[input$simdatadir]]$CX_READ_LENGTH)),
					sort(unique(md[[input$simdatadir]]$CX_READ_LENGTH))),
				checkboxGroupInput("depth", "Read Depth (coverage)",
					sort(unique(md[[input$simdatadir]]$CX_READ_DEPTH)),
					sort(unique(md[[input$simdatadir]]$CX_READ_DEPTH))),
				checkboxGroupInput("fragsize", "Mean Library Fragment Size",
					sort(unique(md[[input$simdatadir]]$CX_READ_FRAGMENT_LENGTH)),
					sort(unique(md[[input$simdatadir]]$CX_READ_FRAGMENT_LENGTH))),
				checkboxGroupInput("aligner", "Aligner",
					PrettyAligner(input$simdatadir),
					"best"),
				checkboxGroupInput("caller", "Software",
					sort(as.character(unique(StripCallerVersion(md[[input$simdatadir]]$CX_CALLER)))),
					sort(as.character(unique(StripCallerVersion(md[[input$simdatadir]]$CX_CALLER)))))
			))
	})
	#####
	# sim plots
	output$simEventSizePlot <- renderPlot({
		cachedsimdata <- RefreshSimData(input, cachedsimdata)
		plotdf <- PrettyFormatSimPlotdf(input, cachedsimdata, cachedsimdata$dfs$callsByEventSize)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf) +
			aes(group=paste(Id, CallSet), x=abs(svLen), y=sens) +
			aes_string(linetype=paste0("as.factor(", input$linetype, ")"), colour=paste0("as.factor(", input$colour, ")")) +
			geom_line(size=0.5) +
			scale_x_svlen +
			# simfacets for the fields not displayed in linetype or colour
			facet_grid(eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$linetype, input$colour))], collapse=" + "))))) +
			labs(title="", y="Sensitivity", x="Event size",
				linetype=names(simfacets)[simfacets==input$linetype],
				colour=names(simfacets)[simfacets==input$colour])
		return(p)
	})
	output$simRocPlot <- renderPlot({
		cachedsimdata <- RefreshSimData(input, cachedsimdata)
		plotdf <- PrettyFormatSimPlotdf(input, cachedsimdata, cachedsimdata$dfs$roc)
			if (nrow(plotdf) == 0) return(NULL)
			p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
				aes(group=paste(Id, CallSet), y=sens, x=fp+1) +
				aes_string(linetype=paste0("as.factor(", input$linetype, ")"), colour=paste0("as.factor(", input$colour, ")")) +
				geom_line() +
				scale_x_log_fp +
				facet_grid(eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$linetype, input$colour))], collapse=" + "))))) +
				labs(title="", y="Sensitivity", x="False Positives",
					linetype=names(simfacets)[simfacets==input$linetype],
					colour=names(simfacets)[simfacets==input$colour])
			return(p)
		})
	output$simbpErrorDistributionPlot <- renderPlot({
		cachedsimdata <- RefreshSimData(input, cachedsimdata)
		plotdf <- PrettyFormatSimPlotdf(input, cachedsimdata, cachedsimdata$dfs$bpErrorDistribution)
		# TODO: incorporate error margin only for truth calls
		p <- ggplot(plotdf %>%
				filter(CX_READ_DEPTH==30, CallSet=="High & Low confidence") %>%
				inner_join(mostSensitiveAligner)) +
			aes(x=bperror, y=rate) +
			geom_bar(stat="identity") +
			facet_grid(caller ~ eventtype)
		return(p)
	})
	#####
	# long read plots
	output$lrPrecRecallPlot <- renderPlot({
		cachedlrdata <- RefreshlrData(input, cachedlrdata)
		plotdf <- cachedlrdata$dfs$roc #PrettyFormatSimPlotdf(input, cachedsimdata, cachedlrdata$dfs$roc)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
						aes(group = paste(Id, CallSet), y = sens, x = fp + 1) +
						geom_line() +
						scale_x_log_fp +
						labs(title = "", y = "Sensitivity", x = "False Positives")
		return(p)
	})
}
