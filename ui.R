
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(
    
    # APP TITLE
    titlePanel("Variant allele fraction"),
    
    # SIDEBAR WITH VARIOUS CONTROLS
    sidebarLayout(
        sidebarPanel(
            
            # NUMBER OF VARIANTS
            numericInput("total", 
                         "Total variants", 500000, 
                         min = 1, max = 100000000, step = 10000),
            
            # HETEROZYGOUS FRACTION
            sliderInput("hetfrac",
                        "Fraction heterozygous",
                        min = 0, max = 1, value = 0.8, step = 0.01, round = 2),
            
            # READ DEPTH
            sliderInput("depth",
                        "Read depth",
                        min = 10, max = 100, value = 60),
            
            # SAMPLE PURITY
            sliderInput("purity",
                        "Sample Purity",
                        min = 0, max = 1, value = 0.7, step = 0.01, round = 2),
            
            # NUMBER OF BINS IN HISTOGRAM
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 100,
                        value = 50)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            
            # RENDER MORE CONTROLS IN SAME DIV AS PLOT
            splitLayout(
                cellWidths = c("25%", "25%", "50%"),
                cellArgs = list(style = "padding: 6px"),
                # MAJOR ALLELE COPY NUMBER STATE
                numericInput("major",
                             "Major allele", 1,
                             min = 0, max = 10, step = 1),
                
                # SERVER GENERATES MINOR ALLELE COPY NUMBER STATE - see server.R
                uiOutput("minor"),
                
                # GERMLINE OR SOMATIC
                selectizeInput("type", "Variant Type", c("Germline", "Somatic"), 
                               options = list(dropdownParent = 'body'))
            ),
            
            # JUST PLOT HISTOGRAM
            plotOutput("distPlot")
            
        )
    )
))