#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("A Perfectly Normal Distribution"),
   
   #Get the parameters for our distribution
   sidebarLayout(
     sidebarPanel(
     
       #mu
       numericInput("mu", h3("Mean"), value = 0),
       
       #sigma^2
       numericInput("sigma2", h3("Variance"), value = 1),
    
       #n
       sliderInput("n", 
                   h3("Number of samples:"),
                   min = 0, 
                   max = 5, 
                   value = 2.5,
                   step = 0.1,
                   pre = '10^')
     ),
      
       # Show a plot of the generated distribution
       mainPanel(plotOutput("distPlot"))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
     
      # generate distribution based on input$mu, input$sigma2, and input$n from ui.R
      x <- rnorm(1:10**input$n, input$mu, input$sigma2)
      
      # draw the histogram of values
      hist(x, col = 4)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

