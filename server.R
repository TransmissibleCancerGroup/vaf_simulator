
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

#' Calculate VAF peak location
#' @param p Tumour purity (float: proportion [0, 1])
#' @param cnt Total chromosome number of segment in tumour (int)
#' @param cnh Total chromosome number of segment in host (int: 0-2)
#' @param nt Number of chromosomes carrying the variant in the tumour (int: nt <= cnt)
#' @param nh Number of chromosomes carrying the variant in the host (int: nh <= cnh)
#' @return VAF peak location (float: proportion [0, 1])
vaf <- function(p, cnt, cnh, nt, nh) {
    stopifnot(p >= 0 & p <= 1)
    stopifnot(nt<=cnt)
    stopifnot(nh<=cnh)
    (nh * (1-p) + nt * p) / (cnt*p + cnh*(1-p))
}

shinyServer(function(input, output) {

    output$minor <- renderUI({

        # GENERATE U.I. FOR MINOR ALLELE COPY NUMBER STATE
        numericInput("minor", "Minor allele", min(input$major, input$minor),
                     min = 0, max = input$major, step = 1)
    })

    output$distPlot <- renderPlot({

        # Process input values
        major.cn <- input$major
        minor.cn <- input$minor
        total.cn <- major.cn + minor.cn
        host.cn <- as.numeric(input$hostcn)

        grid<-expand.grid(host = 0:input$hostcn, tumour = 0:total.cn)
        grid$vaf <- apply(grid, 1, function(r) {
            vaf(input$purity, total.cn, host.cn, r[["tumour"]], r[["host"]])
        })

        if (input$type == "Germline") {
            grid <- grid[grid$host>0 & grid$host <= host.cn, ]
        }
        else {
            grid <- grid[grid$host==0 & grid$tumour > 0, ]
        }

        nvar <- input$total / nrow(grid)
        x <- apply(grid, 1, function(r) {
            rbeta(nvar, rpois(nvar, r[["vaf"]]*input$depth), rpois(nvar, (1-r[["vaf"]])*input$depth))
        })

        # draw the histogram with the specified number of bins
        hist(x, breaks = input$bins, xlim = c(0, 1), freq = FALSE, col = 'darkgray', border = 'white')
    })
})
