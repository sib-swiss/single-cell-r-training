FROM rocker/rstudio:4.5.2

# Set project directory
WORKDIR /project

# Use multiple threads for renv
ENV RENV_INSTALL_THREADS=8
ENV RENV_PATHS_CACHE=/renv_cache

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxt6 \
    libgl1 \
    libcurl4-openssl-dev \
    libglpk-dev \
    libxml2-dev \
    libproj-dev \
    libudunits2-dev \
    libgdal-dev \
    tree \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*

# Create cache directory
RUN mkdir -p /renv_cache

# Copy renv lockfile
COPY renv.lock renv.lock

# Install renv and restore packages
RUN Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")'
RUN Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE)'
