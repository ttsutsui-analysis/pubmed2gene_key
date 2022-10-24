## Installation:

#### If you/your lab already have shinyServer
Install required library  

	apt-get update
	apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
	  	
Install required R packages

	pkgs <- c("rentrez","XML","shinydashboard","stringr","tm","DT","ggplot2","future","future.apply","progressr","pbapply","viridis","dplyr")  
	install.packages(pkgs)
	
Then clone repository  

	cd /srv/shiny-server
	git clone git@github.com:ttsutsui-analysis/pubmed2gene_key.git
	
#### For local machine without shinySever
This requires docker. If you don't have docker, download from [https://docs.docker.com/engine/install/]().  

	docker pull ghcr.io/ttsutsui-analysis/ttsutsui/pub2gene_key:1.0
	docker run --rm -d --name pub2gene_key -p 10020:3838 ghcr.io/ttsutsui-analysis/ttsutsui/pub2gene_key:1.0
You can access index page via [http://localhost:10020]()

## Introduction
### pub2gene_human
 
This app will show genes of which publication enriched in certain search words in pubmed. I believe that this app will be useful when researcher need to know typical gene for certain pubmed search result. So far, this app only adapt human genes only.
Usage introduction written below.
 
1. If you work under proxy server, please fill its URL.
2. Fill search words and select used CPUs and push search button.
![pub2geneFig1](https://user-images.githubusercontent.com/116254113/197516372-fe2d17a1-bae4-4c4a-815c-244aabf1fc51.png)
3. App will search your words in pubmed and obtain all PMIDs. Then calculate overlapping PMIDs for each gene links. The genes which has many overlapping PMIDs with search words up-sift from regression curve. So your interesting genes will show higher standard residual value. You may want to decide cut-off value for standard residual depend on your research interest arbitrarily.
![pub2geneFig2](https://user-images.githubusercontent.com/116254113/197516557-8f925cd6-e507-4767-961a-89db7e5a5381.png)
4. Select your cut-off value for each option and push cut-off button. FYI, since residuals are standarized, 1.96 means top 5%.
![pub2geneFig3](https://user-images.githubusercontent.com/116254113/197518602-2d0b8fca-49c8-4f89-8638-9e6a8cb52c60.png)
5. You can have cut-off result. This result can be downloaded as csv file. In the csv, you can find overlapped PMIDs for each gene with search words.
![pub2geneFig4](https://user-images.githubusercontent.com/116254113/197518874-586c2a46-5a57-4cf9-87f2-b710012e3ec1.png)
  
### pub2keyword_human
