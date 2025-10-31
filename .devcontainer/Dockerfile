# Use a base image with R and RStudio Server (from rocker project)
FROM rocker/rstudio:latest

# Install additions you may need (git, system libs for R packages)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# (Optional) create a project user and set working dir
ENV USER rstudio
WORKDIR /home/rstudio

# Expose the RStudio Server port
EXPOSE 8787

# You can optionally install R packages here or rely on startup script