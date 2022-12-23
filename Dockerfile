FROM rocker/shiny:4.2.1

RUN apt-get update && apt-get install -y libxml2 wget libtime-hires-perl

RUN R -e "options(warn=2); install.packages('rentrez')"
RUN R -e "options(warn=2); install.packages('XML')"
RUN R -e "options(warn=2); install.packages('shinydashboard')"
RUN R -e "options(warn=2); install.packages('DT')"
RUN R -e "options(warn=2); install.packages('ggplot2')"
RUN R -e "options(warn=2); install.packages('future')"
RUN R -e "options(warn=2); install.packages('future.apply')"
RUN R -e "options(warn=2); install.packages('progressr')"
RUN R -e "options(warn=2); install.packages('pbapply')"
RUN R -e "options(warn=2); install.packages('viridis')"
RUN R -e "options(warn=2); install.packages('dplyr')"

RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/versions/18.2.20221128/edirect.tar.gz
RUN tar xzf edirect.tar.gz 
ENV PATH $PATH:/edirect

WORKDIR /edirect
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/versions/18.2.20221128/xtract.Linux.gz 
RUN gunzip -f xtract.Linux.gz  
RUN chmod +x xtract.Linux

RUN rm -rf /srv/shiny-server/*
COPY pub2gene_human /srv/shiny-server/pub2gene_human
COPY pub2keyword_human /srv/shiny-server/pub2keyword_human

