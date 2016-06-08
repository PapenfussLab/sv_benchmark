
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")
source("libplot.R")
library(shiny)

shinyUI(fluidPage(
	titlePanel("Simulated heterzygous simple events"),
	sidebarLayout(
		sidebarPanel(
			selectInput("data", "Data Set", ds),
			#TODO http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
			checkboxInput("smallevents", "Include <= 50bp", value=TRUE),
			uiOutput("simulationControls"),
			hr(),
			selectInput("linetype", "Line Type", facets, "CallSet"),
			selectInput("colour", "Colour", facets, "aligner"),
			width=3
		),
		# Show a plot of the generated distribution
		mainPanel(
			tabsetPanel(
			  tabPanel("Event Size", plotOutput("eventSizePlot", height=1200)),
				tabPanel("ROC",
					plotOutput("rocPlot", height=1200))
			)
		)
	)
))
