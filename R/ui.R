
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")
library(shiny)

shinyUI(fluidPage(
	titlePanel("Heterzygous event simulation"),
	sidebarLayout(
		sidebarPanel(
			selectInput("data", "Data Set", ds),
			conditionalPanel(
				condition="input.data=='rd'",
				selectInput("eventtype", "Event Type",
					withnames(sort(unique(md[["rd"]]$CX_REFERENCE_VCF_VARIANTS)), PrettyVariants(sort(unique(md[["rd"]]$CX_REFERENCE_VCF_VARIANTS)))),
					"hetDEL"),
				checkboxGroupInput("readlength", "Read Length (paired-end)",
					sort(unique(md[["rd"]]$CX_READ_LENGTH)),
					100),
				checkboxGroupInput("depth", "Depth of Coverage",
					sort(unique(md[["rd"]]$CX_READ_DEPTH)),
					60),
				checkboxGroupInput("fragsize", "Mean Library Fragment Size",
					sort(unique(md[["rd"]]$CX_READ_FRAGMENT_LENGTH)),
					300),
				checkboxGroupInput("aligner", "Aligner",
					PrettyAligner("rd"),
					"best"),
				checkboxInput("smallevents", "Include events <= 50bp", value=TRUE)
			),
			#TODO http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale

			#sliderInput("eventsize",
			#	"Log2 event Size:",
			#	min = 0,
			#	max = 16,
			#	value = c(0, 16),
			#	dragRange=TRUE),
			hr(),
			selectInput("linetype", "Line Type", facets, "CallSet"),
			selectInput("colour", "Colour", facets, "aligner"),
			width=3
		),
		# Show a plot of the generated distribution
		mainPanel(
		  plotOutput("mainPlot", height=1000),
			conditionalPanel(
				condition="input.data=='rd'"
			)
		)
	)
))
