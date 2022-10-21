library(shiny)
library(shinydashboard)
library(rentrez)
library(XML)
library(MASS)
library(DT)
library(ggplot2)
library(dplyr)
library(future)

source("functions.R")
options("scipen"=100)

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
  htmlOutput("title1"),
  plotOutput("plot1"),
  DTOutput("table1"),
  p(style = "margin-bottom: 10px;"),
  hr(),
  sliderInput(inputId = "std_res", label="cut-off value for standard residual", min=0,max=10,value=1),
  sliderInput(inputId = "min_ov", label="cut-off value for minimum overlap", min = 1,max = 10,value = 2),
  actionButton(inputId = "Cutoff", label = "Cut-off"),
  br(),
  htmlOutput("title2"),
  plotOutput("plot2"),
  br(),
  DTOutput("table2"),
  downloadButton("df_cut", "Download cut-off result as csv")
)
ui <- dashboardPage(header = header,sidebar =  side,body =  body)

#
server <- function(input,output, session){
  withProgress(message = "Loading",{
    incProgress(0.1)
    gene.all <- readRDS("data/gene_allpub_ids_list.rds")
    incProgress(0.9)})

  reac <- reactiveValues(term=NA, ids=NA, df1=NA, df_cut=NA, proxy=NA, APIkey=NA)
 
  observeEvent(input$APIreg, {
    withProgress(message = "regitering API key", {
    reac$APIkey <- input$APIkey
    for(i in 1:10){Sys.sleep(0.03); incProgress()}
   })
 })
 
  observeEvent(input$proxyreg, {
    withProgress(message = "regitering proxy", {
    reac$proxy <- input$proxy
    for(i in 1:10){Sys.sleep(0.03); incProgress()}
    })
 })
 
  observeEvent(input$Search, {
    if(nchar(reac$APIkey, keepNA = F)==36 & !is.na(reac$APIkey)){
     set_entrez_key(reac$APIkey)
   }
    if(nchar(reac$proxy, keepNA=F)>0 & !is.na(reac$proxy)){
     Sys.setenv(http_proxy=reac$proxy)
     Sys.setenv(https_proxy=reac$proxy)
    }
    
    reac$term <- NA
    output$plot1 <- NULL
    output$table1 <- NULL
    output$plot2 <- NULL
    output$table2 <- NULL
    output$title1 <- NULL
    output$title2 <- NULL
    
    reac$term <- input$query
    if(nchar(reac$term)>0){
      withProgress(message="searching", {
        incProgress(0.3)
        res <- entrez_search(db="pubmed", term=reac$term, retmax=10^6)
        reac$ids <- res$ids
        total <- sapply(gene.all, length)
      })
      
      tf <- ifelse(length(reac$ids)>0,!is.na(reac$ids),FALSE)
      if(tf){
        ids <- reac$ids
        if(input$select==1){
          print("single")
          withProgress(message = "Calculating overlaps",detail = "single mode",{
            count <- count.sapply(gene.all, ids)
          })
        } else {
          withProgressShiny(message = "Calculating overlaps",detail = "multi mode",expr = {
            count <- count.para(all=gene.all, pmid=ids, cores=input$select)
          })
        }
        
        withProgress(message="table formatting",{
          df <- data.frame(symbol=names(gene.all), total_pub=total, overlap_pub=count)
          df <- df[df$total_pub>0, ]
          if(sum(count)>0){
            model <- glm.nb(overlap_pub ~ log10(total_pub), data=df)
            incProgress(0.2)
            std_res <- rstandard(model)
            x <- c(2*4^(0:7))
            y <- predict(model, newdata = data.frame(total_pub=x))
            y <- exp(y)
            
            df.pred <- data.frame(x=x,y=y)
            
            std_res <- round(std_res, 3)
            updateSliderInput(session=session, inputId = "std_res", 
                              min=0,max=max(std_res), value = 1.5, step = 0.01)
            incProgress(0.4)
            df$std_res <- std_res
            df <- df[order(df$std_res, decreasing = T), ]
            reac$df1 <- df
            
            output$title1 <- renderText({print("<font size='+2'><b> Search result </b></font>")})
            output$table1 <- renderDT({datatable(df, filter="top", rownames = F)})
            
            g1 <- ggplot(df, aes(x=total_pub, y=overlap_pub))
            g1 <- g1 + geom_point(aes(color=std_res))
            g1 <- g1 + scale_x_log10() + scale_y_log10() 
            g1 <- g1 + scale_color_gradientn(colours = viridis::turbo(8))
            g1 <- g1 + geom_smooth(data=df.pred, aes(x=x,y=y), method = "lm")
            g1 <- g1 + coord_cartesian(ylim = c(0.9, max(count)*1.1))
            g1 <- g1 + labs(x="total publication for each gene",y="overlapped publication with search result")
            g1 <- g1 + annotate("text", x=1, y=max(df$overlap_pub), label=paste("total number of results=", length(reac$ids)),hjust=0,size=5)
            g1 <- g1 +theme(axis.text=element_text(size=12,family = "sans"),
                            axis.title=element_text(size=14,face="bold"),
                            plot.title = element_text(size=14, face="bold"))
            output$plot1 <- renderPlot({plot(g1)})
            incProgress(0.4)
          } else {
            output$title1 <- renderText({print("<b><font size='+1'> no overlaped publication with gene </font></b>")})
          }
          })
      } else{
        output$title1 <- renderText({print("<b><font size='+2'> no article found </font></b>")})
      }
    }
  })
  
  observeEvent(input$Cutoff,{
      df<- reac$df1
      df1 <- df
      df1$TF <- df$std_res>= input$std_res & df$overlap_pub>= input$min_ov
      
      g2 <- ggplot(df1, aes(x=total_pub, y=overlap_pub))
      g2 <- g2 + geom_point(aes(color=TF)) + scale_color_manual(values=c("gray60","green2"))
      g2 <- g2 + scale_x_log10() + scale_y_log10()
      g2 <- g2 + labs(x="total publication for each gene",y="overlapped publication with search words", color="cut off genes")
      g2 <- g2 +theme(axis.text=element_text(size=12,family = "sans"),
                      axis.title=element_text(size=14,face="bold"),
                      plot.title = element_text(size=14, face="bold"))
      output$plot2 <- renderPlot({plot(g2)})
      
      df_cut <- df[df1$TF, ]
      
      #ov_ids <- find_ov(gene.all, df_cut$symbol, reac$ids)
      #df_cut$overlaped_pmid <- as.character(ov_ids)

      reac$df_cut <- df_cut
      output$title2 <- renderText({print("<font size='+2'><b> Cut-off result </b></font>")})
      output$table2 <- renderDT({datatable(df_cut, filter="top", rownames = F)})
      output$df_cut <- downloadHandler(
        filename = function() {
          paste0(reac$term, " std_residual=",input$std_res,".csv")
        },
        content = function(file) {
          write.csv(reac$df_cut, file, row.names = F)
        }
      )
    })
}      

shinyApp(ui=ui, server = server)