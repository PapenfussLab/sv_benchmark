source("global.R")
library(shiny)
library(ggplot2)
library(dplyr)

cacheddata <- NULL
LoadData <- function(input) {
	if (input$datasettype == "sim") {
		currentdata <- RefreshSimData(input, cacheddata)
	} else {
		currentdata <- RefreshlrData(input, cacheddata)
	}
	cacheddata <<- currentdata
	return(currentdata)
}
RefreshSimData <- function(input) {
  write("RefreshSimData", stderr())
	if (!input$simsmallevents) {
		mineventsize <- 51
	} else {
		mineventsize <- 0
	}
	pd <- LoadPlotData(
		datadir=paste0(dataLocation, "data.", input$simdatadir),
		maxgap=as.numeric(input$simmaxgap),
		ignore.strand=input$simignore.strand,
		sizemargin=input$simsizemargin,
		ignore.duplicates=input$simignore.duplicates,
		ignore.interchromosomal=simoptions$ignore.interchromosomal,
		mineventsize=mineventsize,
		maxeventsize=simoptions$maxeventsize,
		requiredHits=simoptions$requiredHits,
		grtransformName="PrimaryHumanOnly",
		grtransform=simoptions$grtransform[["PrimaryHumanOnly"]],
		truthbedpedir=NULL,
		mintruthbedpescore=NULL,
		eventtypes=NULL,
		existingCache=cacheddata)
	return(pd)
}
RefreshlrData <- function(input, olddata) {
  write("RefreshlrData", stderr())
    pd <- LoadPlotData(
        datadir = paste0(dataLocation, "data.", input$lrdatadir),
        maxgap = input$lrmaxgap,
        ignore.strand = input$lrignore.strand,
        sizemargin = input$lrsizemargin,
        ignore.duplicates = input$lrignore.duplicates,
        ignore.interchromosomal = lroptions$ignore.interchromosomal,
        mineventsize = lroptions$mineventsize,
        maxeventsize = lroptions$maxeventsize,
        requiredHits = input$lrrequiredHits,
        grtransformName = input$lrgrtransformName,
        grtransform = lroptions$grtransform[[input$lrgrtransformName]],
        truthbedpedir = paste0(dataLocation, "input.", input$lrdatadir, "/", lroptions$truthpath[[1]]),
        mintruthbedpescore=input$lrmintruthscore,
				eventtypes=input$lrevents,
        existingCache=cacheddata)
    return(pd)
}
PrettyFormatPlotdf <- function(input, dfs, plotdf) {
	if (input$datasettype == "sim") {
		return(PrettyFormatSimPlotdf(input, dfs, plotdf))
	} else {
		return(PrettyFormatLrPlotdf(input, dfs, plotdf))
	}
}
PrettyFormatLrPlotdf <- function(input, dfs, plotdf) {
	plotdf <- plotdf %>%
		filter(StripCallerVersion(CX_CALLER, FALSE) %in% input$lrcaller)
		# aligner filters
		alignerIdCallSet <- plotdf %>%
				select(Id, CallSet, CX_ALIGNER) %>%
				filter(CX_ALIGNER %in% input$lraligner | (is.na(CX_ALIGNER) & "" %in% input$lraligner)) %>%
				select(Id, CallSet) %>%
				rbind(dfs$mostSensitiveAligner[rep("best" %in% input$lraligner, nrow(dfs$mostSensitiveAligner)),]) %>%
				distinct()
		plotdf <- plotdf %>% inner_join(alignerIdCallSet)
		plotdf <- plotdf %>% mutate(
			caller=StripCallerVersion(CX_CALLER, FALSE),
			aligner=CX_ALIGNER %na% rep("N/A", nrow(plotdf)))
	return(plotdf)
}
PrettyFormatSimPlotdf <- function(input, dfs, plotdf) {
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
				rbind(dfs$mostSensitiveAligner[rep("best" %in% input$simaligner, nrow(dfs$mostSensitiveAligner)),]) %>%
				distinct()
		plotdf <- plotdf %>% inner_join(alignerIdCallSet)
		plotdf <- plotdf %>% mutate(
			caller=StripCallerVersion(CX_CALLER, FALSE),
			aligner=CX_ALIGNER %na% rep("N/A", nrow(plotdf)))
	return(plotdf)
}

doPlot <- function(input, plotdfname, plotfunction, debugLabel=NULL) {
	if (!is.null(debugLabel)) {
		write(debugLabel, stderr())
	}
	progress <- shiny::Progress$new()
	on.exit(progress$close())
	progress$set(message = "Loading data", value = 0)
	current <- RefreshData(input)
	progress$set(message = "Formatting data", value = 0.5)
	plotdf <- PrettyFormatLrPlotdf(input, cachedlrdata, cachedlrdata$dfs[["plotdfname"]])
	# display breakpoint counts instead of breakend counts
	progress$set(message = "Generating plot", value = 0.8)
	p <- plotfunction(plotdf)
}

# server function (must be last in file)
function(input, output, session) {
	#####
	# Bookmarking
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
	  write("simEventSizePlot", stderr())
		progress <- shiny::Progress$new()
		on.exit(progress$close())

		progress$set(message = "Loading data", value = 0)
		currentdata <- RefreshData(input)
		progress$set(message = "Formatting data", value = 0.5)

		currentdata <- LoadData(input)
		plotdf <- PrettyFormatPlotdf(input, currentdata$dfs, currentdata$dfs$callsByEventSize)
		if (nrow(plotdf) == 0) {
		  write("simEventSizePlot: no data!", stderr())
			#browser()
		  return(NULL)
		}
		progress$set(message = "Generating plot", value = 0.8)
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
		return(p)
	})
	output$simRocPlot <- renderPlot({
	  write("simRocPlot", stderr())
		currentdata <- LoadData(input)
		plotdf <- PrettyFormatPlotdf(input, currentdata$dfs, currentdata$dfs$roc)
		p <- plotFacetedRoc(plotdf,
			linetype=paste0("as.factor(", input$simlinetype, ")"),
			colour=paste0("as.factor(", input$simcolour, ")"),
			facet_eval_string=eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + ")))))
		p <- p + labs(linetype=names(simfacets)[simfacets==input$simlinetype],
									colour=names(simfacets)[simfacets==input$simcolour])
	})
	output$simbpErrorDistributionPlot <- renderPlot({
	  write("simbpErrorDistributionPlot", stderr())
		currentdata <- LoadData(input)
		plotdf <- PrettyFormatPlotdf(input, currentdata$dfs, currentdata$dfs$bpErrorDistribution)
		# TODO: incorporate error margin only for truth calls
		p <- ggplot(plotdf %>%
				filter(CX_READ_DEPTH==30, CallSet=="All Calls") %>%
				inner_join(mostSensitiveAligner)) +
			aes(x=bperror, y=rate) +
			geom_bar(stat="identity") +
			facet_grid(caller ~ eventtype)
		write("simbpErrorDistributionPlot complete", stderr())
		return(p)
	})
	#####
	# long read plots
	output$lrPrecRecallPlot <- renderPlot({
		return(doPlot(input, "roc", plotPrecRecall, debugLabel="lrPrecRecallPlot"))
	})
	output$lrRocPlot <- renderPlot({
		return(doPlot(input, "roc", function(plotdf) {
			# display breakpoint counts instead of breakend counts
			plotdf$tp <- plotdf$tp/2
			plotdf$fp <- plotdf$fp/2
			if (nrow(plotdf) == 0) return(NULL)
			p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
				aes(group = paste(Id, CallSet), y = tp, x = fp, colour=caller, linetype=CallSet) +
				geom_line(size=1) +
				coord_cartesian(ylim=c(0, max(plotdf$tp)), xlim=c(0, max(plotdf$tp))) +
				scale_colour_brewer(palette = "Paired") +
				labs(title = "", y = "True Positives", x = "False positives")


			plotroc(plotdf, colour="caller", linetype="CallSet")
		}, debugLabel="lrRocPlot"))
	})
	output$lrPrecRecallRepeatPlot <- renderPlot({
		return(doPlot(input, "rocbyrepeat", plotPrecRecallRepeat, debugLabel="lrPrecRecallRepeatPlot"))
	})
	output$lrRocRepeatPlot <- renderPlot({
		return(doPlot(input, "rocbyrepeat", function(plotdf) {
			plotdf <- PrettyFormatLrPlotdf(input, cachedlrdata, cachedlrdata$dfs$rocbyrepeat)
			# display breakpoint counts instead of breakend counts
			plotdf$tp <- plotdf$tp/2
			plotdf$fp <- plotdf$fp/2
			if (nrow(plotdf) == 0) return(NULL)
			p <- ggplot(plotdf %>% arrange(desc(QUAL))) +
				aes(group = paste(Id, CallSet), y = tp, x = fp, colour=caller, linetype=CallSet) +
				geom_line(size=1) +
				facet_wrap(~ repeatClass) +
				coord_cartesian(ylim=c(0, max(plotdf$tp)), xlim=c(0, max(plotdf$tp))) +
				scale_colour_brewer(palette = "Paired") +
				labs(title = "", y = "True Positives", x = "False positives")
			return(p)
		}, debugLabel="lrRocRepeatPlot"))
	})
}
