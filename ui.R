
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
                         "Total variants", 50000,
                         min = 0, max = 100000000, step = 10000),

            # READ DEPTH
            sliderInput("depth",
                        "Read depth",
                        min = 10, max = 1000, value = 60),
            
            # HOST COPY NUMBER
            radioButtons("hostcn",
                        "Host Copy Number",
                        choices = c(1, 2), selected = 2, inline = TRUE),

            # SAMPLE PURITY
            sliderInput("purity",
                        "Sample Purity",
                        min = 0, max = 1, value = 0.7, step = 0.01, round = 2),

            # FALSE POSITIVES / FALSE NEGATIVES
            splitLayout(
                cellWidths = c("50%", "50%"),
                cellArgs = list(style = "padding: 6px"),
                sliderInput("fpos",
                            "False Positives",
                            min = 0, max = 0.1, value = 0.01, step = 0.01, round = 2),
                sliderInput("fneg",
                            "False Negatives",
                            min = 0, max = 0.1, value = 0.01, step = 0.01, round = 2)
            ),

            # NUMBER OF BINS IN HISTOGRAM
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 100,
                        value = 50),

            # PEAK HEIGHTS
            uiOutput("peak_heights")
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
                selectizeInput("type", "Variant Type", c("Germline", "Somatic", "Both"),
                               options = list(dropdownParent = 'body'))  # So the drop down renders on top of background
            ),

            # JUST PLOT HISTOGRAM
            plotOutput("distPlot", hover = hoverOpts("plot_hover", delay=10)),
            
            # DISPLAY HOVER INFO
            tableOutput("grid")
        )
    )
))
