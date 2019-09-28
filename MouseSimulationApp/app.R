#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(raster)
library(Rcpp)
sourceCpp('Backend.cpp')

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Catch Density"), 
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("iterations", "Iterations of Experiment", min = 1, max = 1000, value = 1),
         sliderInput("ts", "Trap Spacing", min = 0.5, max = 5, value = 1.5),
         sliderInput("delta", "Catch Radius", min = 0.25, max = 2.5, value = 0.5),
         sliderInput("density", "Population Density", min = 0.5, max = 5, value = 1),
         sliderInput("boarder", "Boarder", min = 1, max = 10, value = 3)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         htmlOutput("errorMessage")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   NUM_SQUARES <- 8
  
   output$distPlot <- renderPlot({
     fs = (input$ts*7)+(2*input$boarder)
     np = input$density * fs * fs
     reps <- as.integer(input$iterations)
     # initial matrix n = 1, then sum the others 
     # GenTrapData returns a trap x day counts 
     trapRes <- GenTrapData(trapSpacing = input$ts, catchRadius = input$delta, boarder = input$boarder, nSquares = NUM_SQUARES, trueDensity = input$density, nForages = 4)
     if (reps > 1) {
        for (i in 2:reps) {
           trapRes <- trapRes + GenTrapData(trapSpacing = input$ts, catchRadius = input$delta, boarder = input$boarder, nSquares = NUM_SQUARES, trueDensity = input$density, nForages = 4)
        }
     }
     # Collect results and plot
     trap_mat <- matrix(rowSums(trapRes) / sum(trapRes), nrow = 2*NUM_SQUARES, ncol = 2*NUM_SQUARES)
     trap_rast <- raster(xmn = 0, xmx = 2*NUM_SQUARES, ymn = 0, ymx = 2*NUM_SQUARES, nrows = 2*NUM_SQUARES, ncols = 2*NUM_SQUARES)
     trap_rast[] <- trap_mat
     plot(trap_rast)
   })
   
   output$errorMessage <- renderText({
     if (input$delta > input$ts/2) {
       paste("<font color='red' size='4'><b>", "WARNING: If the catch radius is greater than half the trap spacing a mouse can be exposed to multiple traps. We do not choose a random trap in this case, and so the results may be skewed.", "</b></font>")
     }
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

