FROM rocker/rstudio:4.5.2

WORKDIR /project

# ----------------------------
# renv settings
# ----------------------------
ENV RENV_INSTALL_THREADS=8
ENV RENV_PATHS_CACHE=/renv_cache

# ----------------------------
# Install system dependencies
# ----------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libnlopt-dev \
    libxt6 \
    libgl1 \
    libcurl4-openssl-dev \
    libssl-dev \
    libglpk-dev \
    libxml2-dev \
    libproj-dev \
    libudunits2-dev \
    libgdal-dev \
    libcairo2-dev \
    libxt-dev \
    tree \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /renv_cache

# ----------------------------
# Copy renv lockfile
# ----------------------------
COPY renv.lock renv.lock

# ----------------------------
# Install renv and restore packages
# ----------------------------
RUN Rscript -e 'install.packages("renv", repos="https://packagemanager.posit.co/cran/latest")' \
    && Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE)'