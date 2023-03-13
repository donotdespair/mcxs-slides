
library(shiny)

# Define colors
mcxs1  = "#05386B"
mcxs2  = "#379683"
mcxs3  = "#5CDB95"
mcxs4  = "#8EE4AF"
mcxs5  = "#EDF5E1"

ui <- fluidPage(
  
  titlePanel("Autocorrelation patterns of AR(1) model"),
  plotOutput("acf"),
  hr(),
  fluidRow(
    column(12, align="center",
           h4("Select the values of the parameter"),
    ),
    fluidRow(12),
    fluidRow(
      column(4, offset=4,
             sliderInput("a1",
                         HTML("autoregressive parameter &alpha;<sub>1</sub>"),
                         min = -1,
                         max = 1,
                         step = 0.05,
                         value = 0.5)
      )
    )
  )
)

server <- function(input, output) {
  
  output$acf <- renderPlot({
    acfs     = rep(NA,20) 
    a1       = input$a1
    acfs[1]  = a1
    for (i in 2:20){
      acfs[i] = a1*acfs[i-1]
    }
    pacf     = c(a1,rep(0,19))
    par(mfrow=c(1,2))
    plot(x=1:20, y=acfs, type="h", ylim=c(-1,1), xlab="lags", ylab="acf", col=mcxs1, bty="l",lwd=10)
    abline(h=0)
    plot(x=1:20, y=pacf, type="h", ylim=c(-1,1), xlab="lags", ylab="pacf", col=mcxs2, bty="l",lwd=10)
    abline(h=0)
  })
}

shinyApp(ui = ui, server = server)

