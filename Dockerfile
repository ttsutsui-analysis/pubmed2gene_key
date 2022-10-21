FROM rocker/shiny:4.2.1

RUN apt-get update && apt-get install -y libxml2

RUN R -e "options(warn=2); install.packages('rentrez')"
RUN R -e "options(warn=2); install.packages('XML')"
RUN R -e "options(warn=2); install.packages('shinydashboard')"
RUN R -e "options(warn=2); install.packages('stringr')"
RUN R -e "options(warn=2); install.packages('tm')"
RUN R -e "options(warn=2); install.packages('DT')"
RUN R -e "options(warn=2); install.packages('ggplot2')"
RUN R -e "options(warn=2); install.packages('future')"
RUN R -e "options(warn=2); install.packages('future.apply')"
RUN R -e "options(warn=2); install.packages('progressr')"
RUN R -e "options(warn=2); install.packages('pbapply')"
RUN R -e "options(warn=2); install.packages('viridis')"
RUN R -e "options(warn=2); install.packages('dplyr')"

RUN rm -rf /srv/shiny-server/*
COPY apps/pub2gene_human /srv/shiny-server/pub2gene_human
COPY apps/pub2keyword_human /srv/shiny-server/pub2keyword_human