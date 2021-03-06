library(shinydashboard)
library(shiny)
library(patchwork)
library(ggplot2)
library(cowplot)
#library(shiny)
#library(rbokeh)
library(tidyverse)
source("data/visualisation.R")
#library(Seurat)
options(future.globals.maxSize = Inf)
#memory.limit(100000)
samples <- read_rds("data/samples.rds")
reference <- read_rds("data/Reference.rds")
genes <- unique(c(reference$features.plot))
a <- c("Pct", "pct.exp")
names(a) <- c("Size of ident", "% ident cells with expression")

ui <- dashboardPage(skin = "green",
  dashboardHeader(title = "DevKidCC Kidney Organoid Gene Explorer",
                  titleWidth = 450),
  
  dashboardSidebar(width = 350,
                   
    selectInput('features', 'Genes', genes, multiple = T, selected = c("DAPL1", "LYPD1", "JAG1",
                                                                       "GATA3", "SLC12A1", "HNF4A",
                                                                       "LRP2", "CLDN1", "OLFM3", 
                                                                       "NPHS2", "CENPF")),
    
    selectInput('idents', 'Segments', c("all", "nephron", unique(reference$Component)),
                selected=c("nephron"), multiple = T),
    
    selectInput('samples', 'Samples', samples$sample, 
                selected = "Reference", multiple = T),
    
    radioButtons('show', 'Plot by Gene or Identity', c("gene", "identity"),
                 selected = "gene", inline = T),
    
    radioButtons('size', 'Dot size', a, selected = "Pct", inline=T),
    
    selectInput('columns', 'Columns', 1:10, selected = 1),
    
    radioButtons("scaling", "Scaling", c("log", "scale", "raw"), selected = "log", inline = T),
    
    sliderInput("dot.min", "Min % Cells", min = 0, max = 25, step = 0.5, value = 1),
    
    sliderInput("dotsize", "Size", min = 1, max = 30, step = 1, value = 6),
   
    submitButton(text = "Render Plot", icon = icon("share"))
  ),
  
  dashboardBody(
              fluidRow(
                column(width = 12,
                       plotOutput("distPlot", height = 1200),
                ),
              
    )
  )
)


server <- function(input, output) {
  plotObject <- reactive({
    comp.data <- map_df(input$samples[input$samples!="Reference"],
                        ~read_rds(paste0("data/", .x, ".rds")))
    comp.data <- bind_rows(comp.data, (reference %>% filter(Tier == "DKCC")))
    comp.data <- comp.data %>% filter(features.plot %in% input$features)
  })
  
  output$distPlot <- renderPlot({
    DotPlotCompare(plotObject(), 
                   features = input$features, split.by = input$sample, 
                   classification = "DKCC",
                   show = input$show, 
                   scaling = input$scaling, idents = input$idents, 
                   size = input$size, dot.min = input$dot.min, 
                   dot.scale = input$dotsize,
                   samples = input$samples,
                   columns = input$columns)
  })
}

shinyApp(ui, server)