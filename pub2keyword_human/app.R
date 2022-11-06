library(shiny)
library(shinydashboard)
library(rentrez)
library(XML)
library(stringr)
library(tm)
library(DT)
library(ggplot2)
library(future)

source("functions.R")
options("scipen"=100)

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
  tags$head(tags$style(".shiny-notification {position: fixed; top: calc(50%); left: calc(50%);width: 30%; max-width: 450px; margin-left: auto;margin-right: auto;")),
  htmlOutput("title1"),
  DTOutput("table1"),
  p(style = "margin-bottom: 10px;"),
  hr(),
  htmlOutput("title2"),
  DTOutput("table2"),
  p(style = "margin-bottom: 10px;"),
  hr(),
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

ui <- dashboardPage(header = header,sidebar =  side,body =  body)

server <- function(input, output, session){
  withProgress(message = "Loading",{
    incProgress(0.1)
    bg.key.df <- readRDS("data/gene_allpub_keyword_df_cleaned.rds")
    bg.key.df <- bg.key.df[order(bg.key.df$Freq, decreasing = T), ]
    incProgress(0.9)})
  
  reac <- reactiveValues(term=NA, ids=NA,keywords=NA,pub_df=NA, key_df=NA, proxy=NA, APIkey=NA)
  
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
    
    output$table1 <- NULL
    output$table2 <- NULL
    output$plot1 <- NULL
    output$title1 <- NULL
    output$title2 <- NULL
    
    reac$term <- input$query
    
    if(nchar(reac$term)>0){
      withProgress(message="searching", {
        incProgress(0.2)
        res <- try({entrez_search(db="pubmed",retmax=10^3,term = reac$term)})
        if(length(res)==1){
          ids <- "connection error"
        } else {
          reac$ids <- res$ids
          incProgress(0.4)
        }
      })
      
      tf <- ifelse(length(reac$ids)>0,!is.na(reac$ids),FALSE)
      if(tf){
        ids <- reac$ids
        withProgressShiny(rec.ls <- get_xml(ids), message="retriving articles", detail =paste0("number of articles=", length(ids)))
        
        withProgress(message = "formatting data", {
          keywords <- sapply(rec.ls, xml2key)
          keywords <- unlist(keywords)
          reac$keywords <- keywords
          incProgress(amount = 0.1)
          tl <- sapply(rec.ls, function(x){
            xpathSApply(x, "/PubmedArticleSet/*", get_title)
          })
          tl <- as.character(unlist(tl))
          incProgress(amount = 0.4)
          key <- sapply(rec.ls, function(x){
            xpathSApply(x, "/PubmedArticleSet/*", get_key)
          })
          key <- as.character(unlist(key))
          incProgress(amount=0.4)
          title.key.df <- data.frame(pmid=ids, title=tl, keywords=key)
          reac$pub_df <- title.key.df
          
          dt <- datatable(title.key.df, filter="top", rownames = F)
          output$table1 <- renderDT({dt})
          output$title1 <- renderText({print("<font size='+2'><b> Pubmed search result </b></font>")})
          incProgress(0.1)
        })
        
        if(nrow(reac$pub_df)>0 & length(reac$keywords)>0){
          
          ## keyword enrichment
          withProgress(message="counting obtained keywords",{
            keywords <- reac$keywords
            term.freq <- as.data.frame(table(keywords))
            term.freq <- term.freq[order(term.freq$Freq, decreasing = T), ]
          })
          
          term.freq.cl <- term.freq.clean(term.freq, bg.key.df$keyword.all, input$select)
          res <- fexac_test(term.freq.cl, bg.key.df, input$select)
          reac$key_df <- res
          
          df.tmp <- res
          df.tmp <- df.tmp[order(df.tmp$Freq_in_query, decreasing = T), ]
          df.tmp <- df.tmp[order(df.tmp$BF_p), ]
          df.tmp$FoldEnrich <- round(df.tmp$FoldEnrich, 1)
          df.tmp$BF_p <- round(df.tmp$BF_p,8)
          
          #updateSliderInput(session, "FE", max=max(df.tmp$FoldEnrich))
          #updateSliderInput(session, "Freq", value=freq_cut)
          names(df.tmp) <- c("Keyword","Fold enrichment","Frequency in query", "Bonferroni corrected pvalue")
          
          output$table2 <- renderDT({datatable(df.tmp, filter="top", rownames = F)})
          output$title2 <- renderText({print("<font size='+2'><b>Enrichment analysis on each keywords</b><font>")})
          output$key_df <- downloadHandler(
            filename = function() {
              paste0(reac$term, " enrichment analysis.csv")
            },
            content = function(file) {
              write.csv(df.tmp, file, row.names = F)
            }
          )
          
          withProgress(message= "plotting", {
            res <- reac$key_df
            df <- data.frame(key=res$Keyword, BF_p=-log10(ifelse(res$BF_p==0, min(res$BF_p[res$BF_p!=0]), res$BF_p)), freq=res$Freq_in_query,FE=res$FoldEnrich)
            df <- df[order(df$freq, decreasing = T), ]
            df <- df[order(df$BF_p, decreasing = T), ]
            df$key <- factor(df$key, levels=rev(df$key))
            
            df <- dplyr::filter(df, BF_p >= -log10(input$BF_p), FE >= as.integer(input$FE), freq >= as.integer(input$Freq))
            
            df <- head(df, input$numkey)
            incProgress(0.5)
            
            if(nrow(df)>0){
              p1 <- ggplot(df, aes(x=BF_p, y=key)) + theme_classic()
              p1 <- p1 + geom_segment(aes(x=0, xend=BF_p, y=key, yend=key))
              p1 <- p1 + geom_point(aes(color=FE, size=freq)) + scale_size() + scale_colour_gradient(low="gray80",high = "red2")
              p1 <- p1 + theme(axis.ticks.y = element_blank()) 
              if(input$check){
                p1 <- p1 + scale_x_log10(expand = c(0,0), limits = c(1,max(df$BF_p)*1.04),n.breaks=10)
              } else {
                p1 <- p1 + scale_x_continuous(expand = c(0,0), limits = c(0,max(df$BF_p)*1.04))
              }
              
              p1 <- p1 + labs(x="-log10 bonferroni corrected p-value", y="Keyword", title="Enrichment analysis", color="Fold enrichmen", size="Freq in query")
              p1 <- p1 + theme(axis.text=element_text(size=12,family = "sans"),
                               axis.title=element_text(size=13,face="bold"),
                               plot.title = element_text(size=14, face="bold"))
              incProgress(0.5)
              output$plot1 <- renderPlot({plot(p1)})
            }
          })
        }
        
      } else if(ids=="connection error"){
        output$title <- renderText({print("<font size='+1'> entrez connection error. Please waith and try again a little bit later. </font>")})
      }else{
        output$title1 <- renderText({print("<font size='+2'><b> No article for search words </b></font>")})
      }
    }
  })
  
  observeEvent(input$plot, {
    withProgress(message="re-plotting", {
      res <- reac$key_df
      df <- data.frame(key=res$Keyword, BF_p=-log10(ifelse(res$BF_p==0, min(res$BF_p[res$BF_p!=0]), res$BF_p)), freq=res$Freq_in_query,FE=res$FoldEnrich)
      df <- df[order(df$freq, decreasing = T), ]
      df <- df[order(df$BF_p, decreasing = T), ]
      df$key <- factor(df$key, levels=rev(df$key))
      df <- dplyr::filter(df, BF_p >= -log10(input$BF_p), FE >= as.integer(input$FE), freq >= as.integer(input$Freq))
      df <- head(df, input$numkey)
      incProgress(0.5)
      
      if(nrow(df)>0){
        p1 <- ggplot(df, aes(x=BF_p, y=key)) + theme_classic()
        p1 <- p1 + geom_segment(aes(x=0, xend=BF_p, y=key, yend=key))
        p1 <- p1 + geom_point(aes(color=FE, size=freq)) + scale_size() + scale_colour_gradient(low="gray80",high = "red2")
        p1 <- p1 + theme(axis.ticks.y = element_blank()) 
        if(input$check){
          p1 <- p1 + scale_x_log10(expand = c(0,0), limits = c(1,max(df$BF_p)*1.04),n.breaks=10)
        } else {
          p1 <- p1 + scale_x_continuous(expand = c(0,0), limits = c(0,max(df$BF_p)*1.04))
        }
        p1 <- p1 + labs(x="-log10 bonferroni corrected p-value", y="Keyword", title="Enrichment analysis", color="Fold enrichmen", size="Freq in query")
        p1 <- p1 + theme(axis.text=element_text(size=12,family = "sans"),
                         axis.title=element_text(size=13,face="bold"),
                         plot.title = element_text(size=14, face="bold"))
        incProgress(0.5)
        output$plot1 <- renderPlot({plot(p1)})
      }
      
    })
    
  })
}

shinyApp(ui,server)
#runApp(shinyApp(ui,server), launch.browser = T)