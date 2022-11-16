library(shiny)
library(shinydashboard)
library(DT)
library(future)

header <- dashboardHeader(title = "Pubmed to genes")

side <- dashboardSidebar(
  br(),
  textInput(inputId = "query",label = "Pubmed search",value = "", width = "300px"),
  selectInput(inputId = "select", label="Used number of threads", choices = 1:availableCores(), selected = 1),
  helpText("If possible use threads more than 4, otherwise single core may be faster"),
  actionButton(inputId = "Search", label="Start search"),
  br(),
  hr(),
  textInput(inputId = "APIkey", label="If you have Entrez API key, fill below", value="API key", width = "200px"),
  actionButton(inputId = "APIreg", label= "register API key"),
  br(),
  hr(),
  textInput(inputId = "proxy", label = "If you work under proxy, fill below", width = "200px", value="http://user:pass@proxyhost:proxyport"),
  actionButton(inputId = "proxyreg",label = "register proxy")
  
)

body <- dashboardBody(
  tags$head(tags$style(".shiny-notification {position: fixed; top: calc(50%); left: calc(40%);width: 30%; max-width: 450px; margin-left: auto;margin-right: auto;")),
  htmlOutput("title1"),
  plotOutput("plot1"),
  DTOutput("table1"),
  p(style = "margin-bottom: 10px;"),
  hr(),
  sliderInput(inputId = "std_res", label="Standard residual cut-off", min=0,max=10,value=1),
  sliderInput(inputId = "min_ov", label="Minimum overlap cut-off", min = 1,max = 10,value = 2),
  actionButton(inputId = "Cutoff", label = "Cut-off"),
  br(),
  htmlOutput("title2"),
  plotOutput("plot2"),
  br(),
  DTOutput("table2"),
  downloadButton("df_cut", "Download cut-off result as csv")
)

shinyUI(dashboardPage(header = header,sidebar =  side,body =  body))