
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(
	titlePanel("Heterzygous event simulation"),
	sidebarLayout(
		sidebarPanel(
			selectInput("data", "Data Set", c("Read Depth", "Read Length", "Fragment Size")),
			checkboxGroupInput("eventtype", "Event Types",
				c("Deletion", "Insertion", "Tandem Duplication", "Inversion", "Breakpoint", "Breakpoint at SINE/ALU"),
				c("Deletion", "Insertion", "Tandem Duplication", "Inversion", "Breakpoint", "Breakpoint at SINE/ALU")),
			sliderInput("eventsize",
				"Log2 event Size:",
				min = 0,
				max = 16,
				value = c(0, 16),
				dragRange=TRUE)
		),
		# Show a plot of the generated distribution
		mainPanel(
			plotOutput("mainPlot")
		)
	)
))
