source("global.R")
library(shiny)
library(ggplot2)
library(dplyr)

cachedsimdata <- NULL
cachedlrdata <- NULL
RefreshSimData <- function(input, olddata) {
	if (!input$simsmallevents) {
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
		requiredHits=simoptions$requiredHits,
		grtransform=simoptions$grtransform,
		truthgr=NULL,
		eventtypes=NULL,
		existingCache = olddata)
}
RefreshlrData <- function(input, olddata) {
    pd <- LoadPlotData(
        datadir = paste0(dataLocation, "data.", input$lrdatadir),
        maxgap = lroptions$maxgap,
        ignore.strand = lroptions$ignore.strand,
        sizemargin = lroptions$sizemargin,
        ignore.duplicates = lroptions$ignore.duplicates,
        ignore.interchromosomal = lroptions$ignore.interchromosomal,
        mineventsize = lroptions$mineventsize,
        maxeventsize = lroptions$maxeventsize,
        requiredHits = lroptions$requiredHits,
        grtransform = lroptions$grtransform[[input$lrblacklist]],
        truthgr = LoadLongReadTruthgr(paste0(dataLocation, "input.", input$lrdatadir, "/longread")),
    	eventtypes=input$lrevents,
        existingCache = olddata)
}
PrettyFormatLrPlotdf <- function(input, lrdata, plotdf) {
	plotdf <- plotdf %>%
		filter(StripCallerVersion(CX_CALLER, FALSE) %in% input$lrcaller)
		# aligner filters
		alignerIdCallSet <- plotdf %>%
				select(Id, CallSet, CX_ALIGNER) %>%
				filter(CX_ALIGNER %in% input$lraligner | (is.na(CX_ALIGNER) & "" %in% input$lraligner)) %>%
				select(Id, CallSet) %>%
				rbind(lrdata$dfs$mostSensitiveAligner[rep("best" %in% input$lraligner, nrow(lrdata$dfs$mostSensitiveAligner)),]) %>%
				distinct()
		plotdf <- plotdf %>% inner_join(alignerIdCallSet)
		plotdf <- plotdf %>% mutate(
			caller=StripCallerVersion(CX_CALLER, FALSE),
			aligner=CX_ALIGNER %na% rep("N/A", nrow(plotdf)))
	return(plotdf)
}
PrettyFormatSimPlotdf <- function(input, simdata, plotdf) {
	plotdf <- plotdf %>%
			filter(
				StripCallerVersion(CX_CALLER, FALSE) %in% input$simcaller &
				CallSet %in% input$simcallset &
				CX_REFERENCE_VCF_VARIANTS %in% input$simeventtype &
				CX_READ_LENGTH %in% as.numeric(input$simreadlength) &
				CX_READ_DEPTH %in% as.numeric(input$simdepth) &
				CX_READ_FRAGMENT_LENGTH %in% as.numeric(input$simfragsize))
		# aligner filters
		alignerIdCallSet <- plotdf %>%
				select(Id, CallSet, CX_ALIGNER) %>%
				filter(CX_ALIGNER %in% input$simaligner | (is.na(CX_ALIGNER) & "" %in% input$simaligner)) %>%
				select(Id, CallSet) %>%
				rbind(simdata$dfs$mostSensitiveAligner[rep("best" %in% input$simaligner, nrow(simdata$dfs$mostSensitiveAligner)),]) %>%
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
	#####
	# sim
	output$simControls <- renderUI({
		return(span(
			selectInput("simeventtype", "Event Type",
					withnames(sort(unique(md[[input$simdatadir]]$CX_REFERENCE_VCF_VARIANTS)), PrettyVariants(sort(unique(md[[input$simdatadir]]$CX_REFERENCE_VCF_VARIANTS)))),
					"hetDEL"),
				checkboxGroupInput("simreadlength", "Read Length (paired-end)",
					sort(unique(md[[input$simdatadir]]$CX_READ_LENGTH)),
					sort(unique(md[[input$simdatadir]]$CX_READ_LENGTH))),
				checkboxGroupInput("simdepth", "Read Depth (coverage)",
					sort(unique(md[[input$simdatadir]]$CX_READ_DEPTH)),
					sort(unique(md[[input$simdatadir]]$CX_READ_DEPTH))),
				checkboxGroupInput("simfragsize", "Mean Library Fragment Size",
					sort(unique(md[[input$simdatadir]]$CX_READ_FRAGMENT_LENGTH)),
					sort(unique(md[[input$simdatadir]]$CX_READ_FRAGMENT_LENGTH))),
				checkboxGroupInput("simaligner", "Aligner",
					PrettyAligner(md[[input$simdatadir]]$CX_ALIGNER),
					"best"),
				checkboxGroupInput("simcaller", "Software",
					sort(as.character(unique(StripCallerVersion(md[[input$simdatadir]]$CX_CALLER)))),
					sort(as.character(unique(StripCallerVersion(md[[input$simdatadir]]$CX_CALLER)))))
			))
	})
	output$simEventSizePlot <- renderPlot({
		cachedsimdata <- RefreshSimData(input, cachedsimdata)
		plotdf <- PrettyFormatSimPlotdf(input, cachedsimdata, cachedsimdata$dfs$callsByEventSize)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf) +
			aes(group=paste(Id, CallSet), x=abs(svLen), y=sens) +
			aes_string(linetype=paste0("as.factor(", input$simlinetype, ")"), colour=paste0("as.factor(", input$simcolour, ")")) +
			geom_line(size=0.5) +
			scale_x_svlen +
			# simfacets for the fields not displayed in linetype or colour
			facet_grid(eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + "))))) +
			labs(title="", y="Sensitivity", x="Event size",
				linetype=names(simfacets)[simfacets==input$simlinetype],
				colour=names(simfacets)[simfacets==input$simcolour])
		browser()
		return(p)
	})
	output$simRocPlot <- renderPlot({
		cachedsimdata <- RefreshSimData(input, cachedsimdata)
		plotdf <- PrettyFormatSimPlotdf(input, cachedsimdata, cachedsimdata$dfs$roc)
			if (nrow(plotdf) == 0) return(NULL)
			p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
				aes(group=paste(Id, CallSet), y=sens, x=fp+1) +
				aes_string(linetype=paste0("as.factor(", input$simlinetype, ")"), colour=paste0("as.factor(", input$simcolour, ")")) +
				geom_line() +
				scale_x_log_fp +
				facet_grid(eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + "))))) +
				labs(title="", y="Sensitivity", x="False Positives",
					linetype=names(simfacets)[simfacets==input$simlinetype],
					colour=names(simfacets)[simfacets==input$simcolour])
			browser()
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
		browser()
		return(p)
	})
	#####
	# long read plots
	output$lrPrecRecallPlot <- renderPlot({
		cachedlrdata <- RefreshlrData(input, cachedlrdata)
		plotdf <- PrettyFormatLrPlotdf(input, cachedlrdata, cachedlrdata$dfs$roc)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
						aes(group = paste(Id, CallSet), y = precision, x = tp, colour=caller, linetype=CallSet) +
						geom_line(size=1) +
		        scale_colour_brewer(palette = "Paired") +
						labs(title = "", y = "Precision", x = "Recall (true positive count)")
		return(p)
	})
	output$lrRocPlot <- renderPlot({
		cachedlrdata <- RefreshlrData(input, cachedlrdata)
		plotdf <- PrettyFormatLrPlotdf(input, cachedlrdata, cachedlrdata$dfs$roc)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
						aes(group = paste(Id, CallSet), y = tp, x = fp, colour=caller, linetype=CallSet) +
						geom_line(size=1) +
						coord_cartesian(ylim=c(0, max(plotdf$tp)), xlim=c(0, max(plotdf$tp))) +
		        scale_colour_brewer(palette = "Paired") +
						labs(title = "", y = "True Positives", x = "False positives")
		return(p)
	})
	output$lrPrecRecallRepeatPlot <- renderPlot({
		cachedlrdata <- RefreshlrData(input, cachedlrdata)
		plotdf <- PrettyFormatLrPlotdf(input, cachedlrdata, cachedlrdata$dfs$rocbyrepeat)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
						aes(group = paste(Id, CallSet), y = precision, x = tp, colour=caller, linetype=CallSet) +
						geom_line(size=1) +
						facet_wrap(~ repeatClass, scales="free") +
		        scale_colour_brewer(palette = "Paired") +
						labs(title = "Precision-Recall by RepeatMasker repeat class", y = "Precision", x = "Recall (true positive count)")
		return(p)
	})
	output$lrRocRepeatPlot <- renderPlot({
		cachedlrdata <- RefreshlrData(input, cachedlrdata)
		plotdf <- PrettyFormatLrPlotdf(input, cachedlrdata, cachedlrdata$dfs$rocbyrepeat)
		if (nrow(plotdf) == 0) return(NULL)
		p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
						aes(group = paste(Id, CallSet), y = tp, x = fp, colour=caller, linetype=CallSet) +
						geom_line(size=1) +
						facet_wrap(~ repeatClass) +
						coord_cartesian(ylim=c(0, max(plotdf$tp)), xlim=c(0, max(plotdf$tp))) +
		        scale_colour_brewer(palette = "Paired") +
						labs(title = "", y = "True Positives", x = "False positives")
		return(p)
	})
}
