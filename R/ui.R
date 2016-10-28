
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")
source("libplot.R")
library(shiny)

shinyUI(fluidPage(
	titlePanel("Structural Variantion Benchmark Results"),
	sidebarLayout(
		sidebarPanel(
			selectInput("sample", "Sample", ds),
			checkboxGroupInput("events", "Event types",
				choices=c("Deletion", "Insertion", "Inversion", "Tandem Duplication"),
				selected=c("Deletion", "Insertion", "Inversion", "Tandem Duplication")),

			checkboxGroupInput("events", "Event types",
				choices=c("Deletion", "Insertion", "Inversion", "Tandem Duplication"),
				selected=c("Deletion", "Insertion", "Inversion", "Tandem Duplication")),
			# color = caller
			# line type = aligner
			selectInput("aligner", "Aligner", aligner)
			# could facet on repeat region
			# checkboxInput("repeat", "By repeat annotation", value=FALSE)
		),
		mainPanel(
			# TODO what other plots should I have?
			plotOutput("precisionRecallPlot", height=1200),
			plotOutput("rocPlot", height=1200)
			# event size histogram of TP and FP
		)
	),
	sidebarLayout(
		sidebarPanel(
			selectInput("data", "Data Set", ds),
			#TODO http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
			checkboxInput("smallevents", "Include <= 50bp", value=TRUE),
			uiOutput("simulationControls"),
			hr(),
			selectInput("linetype", "Line Type", facets, "CallSet"),
			selectInput("colour", "Colour", facets, "aligner"),
		),
		mainPanel(
			tabsetPanel(
				tabPanel("Event Size", plotOutput("eventSizePlot", height=1200)),
				tabPanel("ROC", plotOutput("rocPlot", height=1200))
			)
		)
	)
))
