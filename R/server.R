source("global.R")
library(shiny)
library(ggplot2)
library(dplyr)

cachedsimdata <- NULL
cachedlrdata <- NULL
RefreshData <- function(input) {
	if (input$datasettype == "sim") {
		cachedsimdata <<- RefreshSimData(input, cachedsimdata)
		return(cachedsimdata)
	} else {
		cachedlrdata <<- RefreshlrData(input, cachedlrdata)
		return(cachedlrdata)
	}
}
RefreshSimData <- function(input, olddata) {
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
		existingCache=olddata)
	return(pd)
}
RefreshlrData <- function(input, olddata) {
  write("RefreshlrData", stderr())
	truthbedpedir <- NULL
	mintruthbedpescore <- NULL
	requiredHits <- 1
	if (input$lrTruthSet == "Long Reads") {
		truthbedpedir <- paste0(dataLocation, "input.", input$lrdatadir, "/", lroptions$truthpath[[1]])
		mintruthbedpescore <- input$lrmintruthscore
		requiredHits <- input$lrrequiredHits
	}
  pd <- LoadPlotData(
      datadir = paste0(dataLocation, "data.", input$lrdatadir),
      maxgap = input$lrmaxgap,
      ignore.strand = input$lrignore.strand,
      sizemargin = input$lrsizemargin,
      ignore.duplicates = input$lrignore.duplicates,
      ignore.interchromosomal = lroptions$ignore.interchromosomal,
      mineventsize = lroptions$mineventsize,
      maxeventsize = lroptions$maxeventsize,
      requiredHits = requiredHits,
      grtransformName = input$lrgrtransformName,
      grtransform = lroptions$grtransform[[input$lrgrtransformName]],
      truthbedpedir = truthbedpedir,
      mintruthbedpescore = mintruthbedpescore,
			eventtypes=input$lrevents,
      existingCache=olddata)
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
			bind_rows(dfs$mostSensitiveAligner[rep("best" %in% input$lraligner, nrow(dfs$mostSensitiveAligner)),]) %>%
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

doPlot <- function(input, plotdfname, plotfunction, debugLabel=NULL, createProgressBar=!is.null(getDefaultReactiveDomain()), data=RefreshData(input)) {
	if (!is.null(debugLabel)) {
		write(debugLabel, stderr())
	}
	if (createProgressBar) {
		progress <- shiny::Progress$new()
		on.exit(progress$close())
		progress$set(message = "Loading data", value = 0)
	}
	if (createProgressBar) {
		progress$set(message = "Formatting data", value = 0.5)
	}
	plotdf <- PrettyFormatPlotdf(input, data$dfs, data$dfs[[plotdfname]])
	# display breakpoint counts instead of breakend counts
	if (createProgressBar) {
		progress$set(message = "Generating plot", value = 0.8)
	}
	p <- plotfunction(plotdf)
	return(p)
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
	output$mainPanelPlots <- renderUI({
		panels <- list()
		if(input$datasettype == "sim") {
			panels <- c(panels, list(
				tabPanel("Event Size", plotOutput("simEventSizePlot", height = mainPlotHeight)),
				tabPanel("ROC", plotOutput("simRocPlot", height = mainPlotHeight))
			))
		}
		panels <- c(panels, list(
			tabPanel("Precision Recall", plotOutput("lrPrecRecallPlot", height = mainPlotHeight)),
			tabPanel("Precision Recall by repeat", plotOutput("lrPrecRecallRepeatPlot", height = mainPlotHeight)),
			#tabPanel("Precision Recall by homology size", plotOutput("lrPrecRecallIhomlenPlot", height = mainPlotHeight)),
			tabPanel("ROC", plotOutput("lrRocPlot", height = mainPlotHeight)),
			tabPanel("ROC by repeat", plotOutput("lrRocRepeatPlot", height = mainPlotHeight)),
			tabPanel("Positional error", plotOutput("bpErrorDistributionPlot", height = mainPlotHeight))
		))
		return(do.call(tabsetPanel, panels))
	})
	output$simEventSizePlot <- renderPlot({
		if (input$datasettype != "sim") {
			write("skipping simEventSizePlot", stderr())
			return (NULL)
		}
		return(doPlot(input, "callsByEventSize", function(plotdf) {
			p <- plotFacettedSensByEventSize(plotdf,
															 linetype=paste0("as.factor(", input$simlinetype, ")"),
															 colour=paste0("as.factor(", input$simcolour, ")"),
															 facet_eval_string=eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + ")))))
			p <- p + labs(linetype=names(simfacets)[simfacets==input$simlinetype],
										colour=names(simfacets)[simfacets==input$simcolour])
			return(p)
		}, debugLabel="simEventSizePlot"))
	})
	output$simRocPlot <- renderPlot({
		if (input$datasettype != "sim") {
			write("skipping simRocPlot", stderr())
			return (NULL)
		}
		return(doPlot(input, "roc", function(plotdf) {
			p <- plotFacetedRocLogFP(plotdf,
													linetype=paste0("as.factor(", input$simlinetype, ")"),
													colour=paste0("as.factor(", input$simcolour, ")"),
													facet_eval_string=eval(parse(text=paste("caller ~ ", paste(simfacets[!(simfacets %in% c(input$simlinetype, input$simcolour))], collapse=" + ")))))
			p <- p + labs(linetype=names(simfacets)[simfacets==input$simlinetype],
										colour=names(simfacets)[simfacets==input$simcolour])
			return(p)
		}, debugLabel="simRocPlot"))
	})
	output$bpErrorDistributionPlot <- renderPlot({
		return(doPlot(input, "bpErrorDistribution", function(plotdf) {
			if (nrow(plotdf) == 0) {
				write("bpErrorDistributionPlot: no data!", stderr())
				return(NULL)
			}
			plotdf <- plotdf %>%
				mutate(fill=if_else(nominalPosition, "Nominal position", "Incorporating caller\nconfidence interval and\nreported microhomology")) %>%
				filter(CallSet==ALL_CALLS)
			p <- ggplot() +
				aes(x=bperror, y=rate, fill=fill) +
				#geom_bar(stat="identity") +
				geom_bar(stat="identity", data=plotdf %>% filter(nominalPosition==FALSE), alpha=0.5) +
				geom_bar(stat="identity", data=plotdf %>% filter(nominalPosition==TRUE), alpha=0.5) +
				scale_fill_brewer(type="qual", palette="Dark2", name="Error margin") +
				facet_wrap( ~ caller) +
				labs(title="Error in called position", x="Error (base pairs)", y="Portion of true positive calls")
			return(p)
		}, debugLabel="bpErrorDistributionPlot"))
	})
	#####
	# long read plots
	output$lrPrecRecallPlot <- renderPlot({
		return(doPlot(input, "roc", plotPrecRecall, debugLabel="lrPrecRecallPlot"))
	})
	output$lrRocPlot <- renderPlot({
		return(doPlot(input, "roc", plotRocLinear, debugLabel="lrRocPlot"))
	})
	output$lrPrecRecallRepeatPlot <- renderPlot({
		return(doPlot(input, "rocbyrepeat", plotPrecRecallRepeat, debugLabel="lrPrecRecallRepeatPlot"))
	})
	#output$lrPrecRecallIhomlenPlot <- renderPlot({
	#	return(doPlot(input, "rocbyihomlen", plotPrecRecallIhomlen, debugLabel="lrPrecRecallIhomlenPlot"))
	#})
	output$lrRocRepeatPlot <- renderPlot({
		return(doPlot(input, "rocbyrepeat", function(plotdf) {
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
