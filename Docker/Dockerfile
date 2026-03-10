# Dockerfile
FROM rocker/rstudio:4.5.2

WORKDIR /project

# ----------------------------
# Use 8 threads for renv and set renv cache path
# ----------------------------
ENV RENV_INSTALL_THREADS=8
ENV RENV_PATHS_CACHE=/renv_cache

# ----------------------------
# Install system dependencies
# ----------------------------
RUN apt-get update && apt-get install -y \
    libxt6 \
    libgl1 \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libxml2-dev \
    libproj-dev \
    libudunits2-dev \
    libgdal-dev \
    tree \
    libcairo2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /renv_cache

# ----------------------------
# Copy renv lockfile
# ----------------------------
COPY renv.lock renv.lock

# ----------------------------
# Install renv and restore packages in a single RUN command
# Ensures renv is in the system library and works on ARM
# ----------------------------
RUN Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")' \
    && Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE)'
