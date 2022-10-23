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

