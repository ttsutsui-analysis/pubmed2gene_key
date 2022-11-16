library(shiny)
library(shinydashboard)
library(DT)
library(future)

header <- dashboardHeader(
  title = "Keyword enrichment analysis")

side <- dashboardSidebar(
  br(),
  textInput(inputId = "query",label = "Pubmed search",value = "", width = "300px"),
  selectInput(inputId = "select", label="Used number of threads", choices = c(1:availableCores()),selected = 1),
  actionButton(inputId = "Search", label="Start search"),
  br(),
  hr(),
  textInput(inputId = "APIkey", label="If you have Entrez API key, fill below", value="API key", width = "200px"),
  actionButton(inputId = "APIreg", label= "register API key"),
  hr(),
  textInput(inputId = "proxy", label = "If you work under proxy, fill below", width = "200px", value="http://user:pass@proxyhost:proxyport"),
  actionButton(inputId = "proxyreg",label = "register proxy")
  
)

body <- dashboardBody(
  tags$head(tags$style(".shiny-notification {position: fixed; top: calc(50%); left: calc(40%);width: 30%; max-width: 450px; margin-left: auto;margin-right: auto;")),
  htmlOutput("title1"),
  DTOutput("table1"),
  p(style = "margin-bottom: 10px;"),
  hr(),
  htmlOutput("title2"),
  DTOutput("table2"),
  p(style = "margin-bottom: 10px;"),
  hr(),
  htmlOutput("title3"),
  plotOutput("plot1"),
  fluidRow(column(4,sliderInput(inputId = "BF_p", label = "Plot cut-off: Bonferroni corrected p-value", min=0,max=0.2,value=0.01, step=0.001),
                  sliderInput(inputId = "FE", label="Plot cut-off: Minimum fold enrichment in query", min=1,max=100, step=1, value=10)
  ),
  column(4,sliderInput(inputId = "Freq", label="Plot cut-off: Minimum frequency in query", min=1, max=50,step=1, value=3),
         sliderInput(inputId = "numkey", label= "Plot cut-off: Number of displayed keyword", min=10, max=50, step=1, value=30))),
  
  checkboxInput(inputId = "check", label = "x-axis log scale"),
  actionButton(inputId = "plot", label = "Re-draw plot"),
  downloadButton("key_df", "Download enrichment result as csv")
)

shinyUI(dashboardPage(header = header,sidebar =  side,body =  body))

