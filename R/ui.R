
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")
source("libplot.R")

# ui function (must be last in file)
function(request) {
	fluidPage(
		titlePanel("Structural Variant Caller Benchmark"),
		sidebarLayout(
			sidebarPanel(
				selectInput("datasettype", "Data set", choices=c("Genome in a Bottle"="lr", "Simulation"="sim"), selected="lr"),
				conditionalPanel("input.datasettype == 'lr'",
					selectInput("lrdatadir", "Sample", lroptions$datadir),
					selectInput("lrevents", "Event types",
						eventtypes,
						"DEL"),
					#selectInput("lrcallset", "Call Set",
					#	c("High confidence only", "High & Low confidence"),
					#	"High & Low confidence"),
					selectInput("lrblacklist", "Blacklist",
            c("None", "ENCODE DAC blacklist"="DAC", "ENCODE Duke blacklist"="Duke"),
            "DAC"),
					conditionalPanel("input.datasettype == 'hidden'",
					  # hide aligner as an option because we haven't realigned the real data with multiple aligners
  					checkboxGroupInput("lraligner", "Aligner",
  						PrettyAligner((bind_rows(md) %>% distinct(CX_ALIGNER) %>% filter(!is.na(CX_ALIGNER)))$CX_ALIGNER),
  						"best")
					),
					checkboxGroupInput("lrcaller", "Software",
						unique(sort(as.character(StripCallerVersion((bind_rows(md) %>% distinct(CX_CALLER) %>% filter(!is.na(CX_CALLER)))$CX_CALLER)))),
						       unique(sort(as.character(StripCallerVersion((bind_rows(md) %>% distinct(CX_CALLER) %>% filter(!is.na(CX_CALLER)))$CX_CALLER)))))
					# line type = call set
					# color = caller
					#selectInput("aligner", "Aligner", knownaligners),
					# could facet on repeat region
					# checkboxInput("repeat", "By repeat annotation", value=FALSE)
				),
				conditionalPanel("input.datasettype == 'sim'",
					selectInput("simdatadir", "Data Set", simoptions$datadir),
					#TODO http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
					checkboxInput("simsmallevents", "Include <= 50bp", value = FALSE),
					uiOutput("simControls"),
					hr(),
					selectInput("simlinetype", "Line Type", simfacets, "CallSet"),
					selectInput("simcolour", "Colour", simfacets, "aligner"),
					checkboxGroupInput("simcallset", "Call Set",
						c("High confidence only", "High & Low confidence"),
						"High & Low confidence")
				),
				submitButton("refresh", "Refresh")
			),
			mainPanel(
			  conditionalPanel("input.datasettype == 'sim'",
			   tabsetPanel(#id = "mtsim",
			     tabPanel("Event Size", plotOutput("simEventSizePlot", height = 1200)),
			     tabPanel("ROC", plotOutput("simRocPlot", height = 1200))
			   )
			  ),
			  conditionalPanel("input.datasettype != 'sim'",
			   tabsetPanel(#id = "mtsim",
			     tabPanel("Precision Recall", plotOutput("lrPrecRecallPlot", height = 1200)),
			     tabPanel("ROC", plotOutput("lrRocPlot", height = 1200)),
			     tabPanel("Precision Recall by repeat", plotOutput("lrPrecRecallRepeatPlot", height = 1200)),
			     tabPanel("ROC by repeat", plotOutput("lrRocRepeatPlot", height = 1200))
			   )
			  )
			  #conditionalPanel("input.datasettype == 'sim'",
					# includeMarkdown("explaination.md")
					# downloadButton('downloadPlot', 'Download Plot')
					# http://stackoverflow.com/questions/14810409/save-plots-made-in-a-shiny-app
			)
		)
	)
}
