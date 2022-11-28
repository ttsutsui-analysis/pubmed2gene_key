library(shiny)
library(rentrez)
library(XML)
library(MASS)
library(DT)
library(ggplot2)
library(dplyr)
library(future)

source("functions.R")
options("scipen"=100)

shinyServer(function(input,output, session){
  withProgress(message = "Loading",{
    incProgress(0.1)
    gene.all <- readRDS("data/gene_allpub_ids_list.rds")
    incProgress(0.9)})
  
  reac <- reactiveValues(term=NA, ids=NA, df1=NA, df_cut=NA, proxy=NA, APIkey=NA, std_res=NA)
  
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
    
    reac$ids <- NA
    reac$term <- NA
    reac$std_res <- NA
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
        res <- try({entrez_search(db="pubmed", term=reac$term, retmax=10^6)})
        if(length(res)==1){
          reac$ids <- "connection error"
        } else {
          if(length(res$ids)>0){
            reac$ids <- res$ids
          } else {
            reac$ids <- NA
          }
          total <- sapply(gene.all, length)
        }
      })
      
      tf <- reac$ids!="connection error" & !is.na(reac$ids)
      if(tf[1]){
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
            std_res <- round(std_res, 3)
            df$std_res <- std_res
            df <- df[order(df$std_res, decreasing = T), ]
            reac$df1 <- df
            updateSliderInput(session=session, inputId = "std_res", 
                              min=0,max=max(std_res), value = 1.5, step = 0.01)
            incProgress(0.4)
            
            x <- c(2*4^(0:7))
            y <- predict(model, newdata = data.frame(total_pub=x))
            y <- exp(y)
            df.pred <- data.frame(x=x,y=y)
            
            g1 <- ggplot(reac$df1, aes(x=total_pub, y=overlap_pub))
            g1 <- g1 + geom_point(aes(color=std_res))
            g1 <- g1 + scale_x_log10() + scale_y_log10() 
            g1 <- g1 + scale_color_gradientn(colours = viridis::turbo(8))
            g1 <- g1 + geom_smooth(data=df.pred, aes(x=x,y=y), method = "lm")
            g1 <- g1 + coord_cartesian(ylim = c(0.9, max(count)*1.1))
            g1 <- g1 + labs(x="total publication for each gene",y="overlapped publication with search words")
            g1 <- g1 + ggplot2::annotate("text", x=1,y=max(reac$df1$overlap_pub),label=paste0("total number of articles \nfor search words=", length(reac$ids)),hjust=0, fontface =2)
            g1 <- g1 +theme(axis.text=element_text(size=12,family = "sans"),
                            axis.title=element_text(size=14,face="bold"),
                            plot.title = element_text(size=14, face="bold"))
            output$plot1 <- renderPlot({plot(g1)})
            output$title1 <- renderText({print("<font size='+2'><b> Search result </b></font>")})
            names(df) <- c("SYMBOL", "total number of pmid for each gene", "number of overlapped pmid", "standard residual")
            output$table1 <- renderDT({datatable(df, filter="top", rownames = F)})
            
          } else {
            output$title1 <- renderText({print("<b><font size='+1'> no overlaped publication with gene </font></b>")})
          }
        })
      } else if(is.na(reac$ids)) {
        output$title1 <- renderText({print("<b><font size='+2'> no article found </font></b>")})
      } else if(reac$ids=="connection error"){
        output$title1 <- renderText({print("<font size='+1'> Entrez connection error.<br>Please wait and try again a little bit later. </font>")})
      }
    }
  })
  
  observeEvent(input$Cutoff,{
    reac$std_res <- input$std_res
    df<- reac$df1
    df1 <- df
    df1$TF <- df$std_res>= input$std_res & df$overlap_pub>= input$min_ov
    
    g2 <- ggplot(df1, aes(x=total_pub, y=overlap_pub))
    g2 <- g2 + geom_point(aes(color=TF)) + scale_color_manual(values=c("gray60","green2"))
    g2 <- g2 + scale_x_log10() + scale_y_log10()
    g2 <- g2 + labs(x="total publication for each gene",y="overlapped publication with search words", color="cut-off genes")
    g2 <- g2 +theme(axis.text=element_text(size=12,family = "sans"),
                    axis.title=element_text(size=14,face="bold"),
                    plot.title = element_text(size=14, face="bold"))
    output$plot2 <- renderPlot({plot(g2)})
    
    df_cut <- df[df1$TF, ]
    
    reac$df_cut <- df_cut
    output$title2 <- renderText({print("<font size='+2'><b> Cut-off result </b></font>")})
    names(df_cut) <- c("SYMBOL", "total number of pmid for each gene", "number of overlapped pmid", "standard residual")
    output$table2 <- renderDT({datatable(df_cut, filter="top", rownames = F)})
    
    output$df_cut <- downloadHandler(
      filename = function() {
        paste0(reac$term, " std_residual=",reac$std_res,".csv")
      },
      content = function(file) {
        withProgress(message = "making csv",{
          df <- reac$df_cut
          symbs <- reac$df_cut$symbol
          df$pmids <- pick.pmid(gene.all, reac$ids, symbs)
          names(df) <- c("SYMBOL", "total number of pmid for each gene", "number of overlapped pmid", "standard residual", "overlapped pmids")
        })
        write.csv(df, file, row.names = F)
      }
    )
  })
})