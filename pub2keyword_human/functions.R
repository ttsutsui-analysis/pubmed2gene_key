library(XML)
library(shiny)
library(pbapply)
library(progressr)
library(future.apply)

# get keyword
get_xml <- function(ids){
  n <- ceiling(length(ids)/200)
  rec.ls <- list()
  keywords <- NA
  p <- progressr::progressor(steps=n)
  for(i in 1:n){
    num <- seq(from=1+200*(i-1), to=200+200*(i-1), by=1)
    ids.i <- ids[num]
    ids.i <- na.omit(ids.i)
    upload <- entrez_post(db="pubmed", id=ids.i)
    rec <- try({entrez_fetch(db="pubmed", rettype = "xml", parsed=TRUE, web_history = upload)})
    if(length(rec)!=1){
      rec.ls <- c(rec.ls, rec)
    } else {
      rec <- entrez_fetch(db="pubmed", rettype = "xml", parsed=TRUE, web_history = upload)
      rec.ls <- c(rec.ls, rec)
    }
    p()
  }
  return(rec.ls)
}

# keyword cleaning tool
xml2key <- function(xml){
  keywords <- xpathSApply(xml, "//KeywordList/Keyword", xmlValue)
  
  keywords <- as.character(unlist(keywords))
  keywords <- strsplit(keywords, split=";")
  keywords <- unlist(keywords)
  keywords <- na.omit(keywords)
  keywords <- gsub("^\\s", "", keywords)
  keywords <- gsub("\\s$", "", keywords)
  keywords <- str_squish(keywords)
  keywords <- removePunctuation(keywords)
  keywords <- gsub("cells", "cell", keywords)
  keywords <- gsub("Cells","Cell", keywords)
  keywords <- gsub("cytokines","cytokine", keywords)
  keywords <- gsub("Cytokines","Cytokine", keywords)
  keywords <- gsub("inhibitors", "inhibitor", keywords)
  keywords <- gsub("Inhibitors", "Inhibitor", keywords)
  return(keywords)
}


# get article title and keyword by article
get_xmlValue<- function(x, node){
  a <- xpathSApply(x, node, xmlValue)
  if(length(a) == 0) {a <- NA} else {a}
}

get_title <- function(x){
  get_xmlValue(x, ".//Article/ArticleTitle")
}

get_key <- function(x){
  key <- get_xmlValue(x, ".//KeywordList/Keyword")
  paste0(key,collapse = ", ")
}

# aggregate uppercase or lowercase into the same
term.freq.clean <- function(df, bg, core){
  withProgress(message = "keyword formatting",{
    core <- as.integer(core)
    df <- df[order(df[,2], decreasing = T), ]
    df2 <- data.frame(term=bg,Freq=NA)
    df <- rbind(df2, setNames(df,names(df2)))
    incProgress()
    terms <- as.character(df[,1])
    # get term stems
    term.stem <- toupper(terms)
    term.stem <- gsub("-", "", term.stem)
    term.stem <- gsub("‐", "", term.stem)
    term.stem <- gsub("–", "", term.stem)
    term.stem <- gsub("\\s","", term.stem)
    incProgress()
    term.stem.dup <- duplicated(term.stem)
    need.to.agg <- term.stem[term.stem.dup]
    need.to.agg <- need.to.agg[!duplicated(need.to.agg)]
    incProgress()
    tmp <- df[!term.stem %in% need.to.agg, ]
    tmp <- na.omit(tmp)
  })
    
  if(core==1){
    withProgress(message="keyword cleansing",detail = "single mode",{
      tmp2 <- pbsapply(need.to.agg, function(x){
        term.i <- x
        df.i <- df[term.stem==term.i, ]
        incProgress(1/length(need.to.agg))
        c(keyword.all=as.character(df.i[1,1]), Freq=sum(df.i$Freq, na.rm = T))
      })
    })
  } else {
    withProgressShiny(message = "keyword cleansing", detail="multi mode",{
      plan(multisession, workers=core)
      p <- progressor(steps = 100)
      n <- floor(length(need.to.agg)/100)
      tmp2 <- future_sapply(seq_along(need.to.agg), function(x){
        term.i <- need.to.agg[x]
        df.i <- df[term.stem==term.i, ]
        if(x %% n ==0)p()
        c(keyword.all=as.character(df.i[1,1]), Freq=sum(df.i$Freq, na.rm = T))
      })
    })
  }

  tmp2 <- data.frame(V1=tmp2[1,], V2=as.integer(tmp2[2,]))
  
  df.ok <- rbind(tmp, setNames(tmp2, names(tmp)))
  
  df.ok <- df.ok[order(df.ok$Freq, decreasing = T), ]
  #df.ok <- df.ok[df.ok$Freq>0, ]
  return(df.ok)
}

# enrichment analysis
fexac_test <- function(key.df, bg.df,core){
  core <- as.integer(core)
  if(core==1){
    withProgress(message="statistical testing", detail = "single mode",{
      res <- pbapply(as.matrix(key.df),1,function(x){
        term <- as.character(x[1])
        first <- c(as.integer(x[2]), sum(key.df$Freq) - as.integer(x[2]))
        if(term %in% bg.df$keyword.all){
          freq.tmp <- bg.df$Freq[bg.df$keyword.all==term]
          second <- c(freq.tmp, sum(bg.df$Freq) - freq.tmp)
        } else {
          second <- c(1, sum(bg.df$Freq)+1)
        }
        p <- fisher.test(rbind(first, second))
        FE <- (first[1]/first[2])/(second[1]/second[2])
        incProgress(1/nrow(key.df))
        c(p$p.value, FE,first[1],second[1],term)
      })
    })
  }else{
    withProgressShiny(message="statistical testing", detail = "multi mode",{
      p <- progressor(steps=100)
      n <- floor(nrow(key.df)/100)
      res <- future_sapply(1:nrow(key.df),function(y){
        x <- key.df[y, ]
        term <- as.character(x[1,1])
        first <- c(as.integer(x[2]), sum(key.df$Freq) - as.integer(x[2]))
        if(term %in% bg.df$keyword.all){
          freq.tmp <- bg.df$Freq[bg.df$keyword.all==term]
          second <- c(freq.tmp, sum(bg.df$Freq) - freq.tmp)
        } else {
          second <- c(1, sum(bg.df$Freq)+1)
        }
        ps <- fisher.test(rbind(first, second))
        FE <- (first[1]/first[2])/(second[1]/second[2])
        
        if(y %% n ==0)p()
        c(ps$p.value, FE, first[1], second[1], term)
      })
      plan(sequential)
    })
  }

  y <- data.frame(Keyword=as.character(res[5,]), FoldEnrich=as.numeric(res[2, ]), 
                  Freq_in_query=as.integer(res[3,]),
                  BF_p=as.numeric(p.adjust(res[1,],method = "bonferroni")))
  row.names(y) <- NULL
  return(y)
}
