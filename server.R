
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

get_grid <- function(input) {
    # Process input values
    major.cn <- input$major
    minor.cn <- input$minor
    total.cn <- major.cn + minor.cn
    host.cn <- as.numeric(input$hostcn)

    grid<-expand.grid(host = 0:host.cn, tumour = 0:total.cn)
    grid$vaf <- apply(grid, 1, function(r) {
        vaf(input$purity, total.cn, host.cn, r[["tumour"]], r[["host"]])
    })
    
    if (input$type == "Germline") {
        grid <- grid[grid$host>0 & grid$host <= host.cn, ]
    }
    else if (input$type == "Somatic") {
        grid <- grid[grid$host==0 & grid$tumour > 0, ]
    }
    grid
}

minor.selected <- 1

shinyServer(function(input, output) {

    observe({
        minor.selected <<- input$minor
    })
    
    output$minor <- renderUI({

        # GENERATE U.I. FOR MINOR ALLELE COPY NUMBER STATE
        numericInput("minor", "Minor allele", min(input$major, minor.selected),
                     min = 0, max = input$major, step = 1)
    })
    
    output$peak_heights <- renderUI({
        req(input$major, input$minor, input$purity)
        # GENERATE U.I. FOR PEAK HEIGHTS
        grid <- grid()
        lapply(1:nrow(grid), function(i) {
            sliderInput(inputId = paste0("Peak", i),
                        label = paste("Peak", sprintf("host%d tumour%d vaf=%.2f", grid[i, 1], grid[i, 2], grid[i, 3])),
                        min = 0, max = 1, value = 0.5, step = 0.01)
        })
    })
    
    grid <- reactive({
        get_grid(input)
    })

    output$distPlot <- renderPlot({
        req(input$major, input$minor, input$purity, input$total)
        grid <- grid()
    
        nvar <- input$total / nrow(grid)
        
        heights <- vector("numeric", nrow(grid))
        
        for (i in 1:nrow(grid)) {
            heights[i] <- input[[paste0("Peak", i)]]
        }
        
        heights <- heights / sum(heights)
        grid$heights <- heights

        x <- unlist(apply(grid, 1, function(r) {
            nv <- round(nvar * r[["heights"]], 0)
            fnr <- input$fneg
            fpr <- input$fpos
            rbeta(nv, rpois(nv, (fpr + r[["vaf"]])*input$depth), rpois(nv, (fnr + 1 - r[["vaf"]])*input$depth))
        }))
        
        # draw the histogram with the specified number of bins
        hist(x, breaks = input$bins, xlim = c(0, 1), freq = TRUE, col = 'darkgray', border = 'white')
        segments(grid$vaf, 0, grid$vaf, -1, lwd = 4, col = "slateblue3")
    })
    
    output$grid <- renderTable({
        # x_str <- function(e) {
        #     if(is.null(e)) return("out of bounds\n")
        #     paste0(round(e$x, 2), "\n")
        # }
        # 
        # paste0("cursor position=", x_str(input$plot_hover))
        grid <- grid()
        heights <- vector("numeric", nrow(grid))
        
        for (i in 1:nrow(grid)) {
            heights[i] <- input[[paste0("Peak", i)]]
        }
        
        heights <- heights / sum(heights)
        grid$heights <- heights
        grid
    })
})

vafhist <- function(nvar, vafs, depth, heights) {
    stopifnot(length(vafs) == length(heights))
    unlist(sapply(1:length(vafs), function(i) {
        nv <- round(nvar*heights[i], 0)
        rbeta(nv, rpois(nv, vafs[i] * depth), rpois(nv, max(0.001, (1-vafs[i])) * depth))
    }))
}