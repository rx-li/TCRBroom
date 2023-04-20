FROM rocker/r-ver:4.2.0
RUN apt-get update && apt-get install -y  jags && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("gpplots",upgrade="never", version = "NA")'
RUN Rscript -e 'remotes::install_version("RColorBrewer",upgrade="never", version = "1.1-3")'
RUN Rscript -e 'remotes::install_version("R2jags",upgrade="never", version = "0.7-1")'
RUN Rscript -e 'remotes::install_version("ComplexHeatmap",upgrade="never", version = "2.12.1")'
RUN R -e 'remotes::install_github("rx-li\/TCRBroom")'
