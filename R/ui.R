
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")
source("libplot.R")

# ui function (must be last in file)
function(request) {
	fluidPage( #TODO try navbarPage
		titlePanel("Structural Variant Caller Benchmark"),
		sidebarLayout(
			sidebarPanel(
				selectInput("lrdatadir", "Sample", lroptions$datadir),
				checkboxGroupInput("events", "Event types",
					choices = c("Deletion", "Insertion", "Inversion", "Tandem Duplication"),
					selected = c("Deletion", "Insertion", "Inversion", "Tandem Duplication")),
				checkboxGroupInput("lrcallset", "Call Set",
					c("High confidence only", "High & Low confidence"),
					"High & Low confidence"),
				# line type = call set
				# color = caller
				#selectInput("aligner", "Aligner", knownaligners),
				# could facet on repeat region
				# checkboxInput("repeat", "By repeat annotation", value=FALSE)
				selectInput("simdatadir", "Data Set", simoptions$datadir),
				#TODO http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
				checkboxInput("simsmallevents", "Include <= 50bp", value = TRUE),
				uiOutput("simControls"),
				hr(),
				selectInput("simlinetype", "Line Type", simfacets, "CallSet"),
				selectInput("simcolour", "Colour", simfacets, "aligner"),
				checkboxGroupInput("simcallset", "Call Set",
					c("High confidence only", "High & Low confidence"),
					"High & Low confidence")
			),
			mainPanel(
				tabsetPanel(id = "mt",
					tabPanel("Precision Recall", plotOutput("lrPrecRecallPlot", height = 1200)),
					tabPanel("ROC", plotOutput("lrocPlot", height = 1200)),
					tabPanel("Detected event sizes", plotOutput("lrEventSizeHistogram", height = 1200)),
					tabPanel("Event Size", plotOutput("simEventSizePlot", height = 1200)),
					tabPanel("ROC", plotOutput("simRocPlot", height = 1200))
					)
				# includeMarkdown("explaination.md")
				# downloadButton('downloadPlot', 'Download Plot')
				# http://stackoverflow.com/questions/14810409/save-plots-made-in-a-shiny-app
			)
		)
	)
}
