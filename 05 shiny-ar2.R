
library(shiny)

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

ui <- fluidPage(
  
  titlePanel("Autocorrelation patterns of AR(2) model"),
  plotOutput("acf"),
  hr(),
  fluidRow(
    column(12, align="center",
           h4("Select the values of the parameters"),
    ),
    fluidRow(12),
    fluidRow(
      column(4, offset=1,
             sliderInput("a1",
                         HTML("autoregressive parameter &alpha;<sub>1</sub>"),
                         min = -1,
                         max = 1,
                         step = 0.05,
                         value = 0.5)
      ),
      column(4, offset=2,
             sliderInput("a2",
                         HTML("autoregressive parameter &alpha;<sub>2</sub>"),
                         min = -1,
                         max = 1,
                         step = 0.05,
                         value = 0.2,
                         uiOutput("a2"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$a1,  {
    updateSliderInput(session = session, inputId = "a2", max = 1 - input$a1)
  })
  observeEvent(input$a2,  {
    updateSliderInput(session = session, inputId = "a1", max = 1 - input$a2)
  })
  
  output$acf <- renderPlot({
    acfs     = rep(NA,20) 
    a1       = input$a1
    a2       = input$a2
    acfs[1]  = (a1)/(1-a2)
    acfs[2]  = (a1^2 - a2^2 + a2)/(1-a2)
    for (i in 3:20){
      acfs[i] = a1*acfs[i-1] + a2*acfs[i-2]
    }
    pacf     = c(a1,a2,rep(0,18))
    par(mfrow=c(1,2))
    plot(x=1:20, y=acfs, type="h", ylim=c(-1,1), xlab="lags", ylab="acf", col=mcxs1, bty="l",lwd=10)
    abline(h=0)
    plot(x=1:20, y=pacf, type="h", ylim=c(-1,1), xlab="lags", ylab="pacf", col=mcxs2, bty="l",lwd=10)
    abline(h=0)
  })
}

shinyApp(ui = ui, server = server)

