
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source("global.R")
source("libplot.R")
mainPlotHeight <- 800
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
					# DAC blacklist was given to some callers so cannot be used generally
					#selectInput("lrblacklist", "Blacklist",
          #  c("None", "ENCODE DAC blacklist"="DAC", "ENCODE Duke blacklist"="Duke"),
          #  "DAC"),
					conditionalPanel("input.datasettype == 'hidden'",
					  # hide aligner as an option because we haven't realigned the real data with multiple aligners
  					checkboxGroupInput("lraligner", "Aligner",
  						PrettyAligner((bind_rows(md) %>% distinct(CX_ALIGNER) %>% filter(!is.na(CX_ALIGNER)))$CX_ALIGNER),
  						"best")
					),
					checkboxGroupInput("lrcaller", "Software",
						unique(sort(as.character(StripCallerVersion((bind_rows(md) %>% distinct(CX_CALLER) %>% filter(!is.na(CX_CALLER)))$CX_CALLER)))),
						       unique(sort(as.character(StripCallerVersion((bind_rows(md) %>% distinct(CX_CALLER) %>% filter(!is.na(CX_CALLER)))$CX_CALLER))))),
					# line type = call set
					# color = caller
					#selectInput("aligner", "Aligner", knownaligners),
					# could facet on repeat region
					# checkboxInput("repeat", "By repeat annotation", value=FALSE)
					hr(),
					helpText("The following inputs control the matching of variant calls against the truth. To reduce bias in the truth set as much as possible, the truth is determined by counting the number of long reads supporting each variant call."),
					selectInput("lrmaxgap", "Maximum break-end error (bp)", lroptions$maxgap, selected=200),
					bsTooltip("lrmaxgap", "Variant calls in which the start or end of the event differs by more than this value are considered unmatched"),
					selectInput("lrsizemargin", "Maximum error in event size", lroptions$sizemargin, selected=0.25),
					bsTooltip("lrsizemargin", "Variant calls in which the event size is greater or smaller than this margin relative to the truth. The default allows up to a 25% difference in event sizes."),
					checkboxInput("lrignore.duplicates", "Ignore duplicate variant calls", lroptions$ignore.duplicates, selected=TRUE),
					bsTooltip("lrignore.duplicates", "Controls whether only the highest scoring variant call matching a given truth call is considered true. The default considers lower scoring duplicate calls as false positives."),
					checkboxInput("lrignore.strand", "Ignore break-end orientation", lroptions$ignore.strand, selected=TRUE),
					bsTooltip("lrignore.strand", "By default, the break-end orientation is ignored since BreakDancer does reports the position by not break-end orientation for variant calls."),
					selectInput("lrgrtransformName", "Blacklist", names(lroptions$grtransform), selected="DAC"),
					bsTooltip("lrgrtransformName", "ENCODE blacklist specifying regions of the genome to ignore. Variant calls with either break-end in a blacklisted region are ignored."),
					selectInput("lrmintruthscore", "Long read MAPQ", lroptions$mintruthscore, selected=0),
					bsTooltip("lrmintruthscore", "Minimum MAPQ score of a long read alignment to be considered as truth."),
					selectInput("lrrequiredHits", "Long read hits", lroptions$requiredHits, selected=5),
					bsTooltip("lrrequiredHits", "Number of long reads required to support the variant to be considered a true positive call.")
				),
				conditionalPanel("input.datasettype == 'sim'",
					selectInput("simdatadir", "Data Set", simoptions$datadir),
					#TODO http://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
					selectInput("simlinetype", "Line Type", simfacets, "CallSet"),
					selectInput("simcolour", "Colour", simfacets, "aligner"),
					checkboxGroupInput("simcallset", "Call Set",
														 c("High confidence only", "High & Low confidence"),
														 "High & Low confidence"),
					checkboxInput("simsmallevents", "Include <= 50bp", value = FALSE),
					uiOutput("simControls"),
					hr(),
					helpText("The following inputs control the matching of variant calls against the truth variants."),
					selectInput("simmaxgap", "Maximum break-end error (bp)", simoptions$maxgap, selected=200),
					bsTooltip("simmaxgap", "Variant calls in which the start or end of the event differs by more than this value are considered unmatched"),
					selectInput("simsizemargin", "Maximum error in event size", simoptions$sizemargin, selected=0.25),
					bsTooltip("simsizemargin", "Variant calls in which the event size is greater or smaller than this margin relative to the truth. The default allows up to a 25% difference in event sizes."),
					checkboxInput("simignore.duplicates", "Ignore duplicate variant calls", simoptions$ignore.duplicates, selected=TRUE),
					bsTooltip("simignore.duplicates", "Controls whether only the highest scoring variant call matching a given truth call is considered true. The default considers lower scoring duplicate calls as false positives."),
					checkboxInput("simignore.strand", "Ignore break-end orientation", simoptions$ignore.strand, selected=TRUE),
					bsTooltip("simignore.strand", "By default, the break-end orientation is ignored since BreakDancer does reports the position by not break-end orientation for variant calls.")
				)
				#,submitButton("refresh", "Refresh")
			),
			mainPanel(
			  conditionalPanel("input.datasettype == 'sim'",
			   tabsetPanel(#id = "mtsim",
			     tabPanel("Event Size", plotOutput("simEventSizePlot", height = mainPlotHeight)),
			     tabPanel("ROC", plotOutput("simRocPlot", height = mainPlotHeight))
			   )
			  ),
			  conditionalPanel("input.datasettype != 'sim'",
			   tabsetPanel(#id = "mtlr",
			     tabPanel("Precision Recall", plotOutput("lrPrecRecallPlot", height = mainPlotHeight)),
			     tabPanel("ROC", plotOutput("lrRocPlot", height = mainPlotHeight)),
			     tabPanel("Precision Recall by repeat", plotOutput("lrPrecRecallRepeatPlot", height = mainPlotHeight)),
			     tabPanel("ROC by repeat", plotOutput("lrRocRepeatPlot", height = mainPlotHeight))
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
