FROM rocker/rstudio:4.6

WORKDIR /project

# ----------------------------
# renv and system settings
# ----------------------------
ENV RENV_INSTALL_THREADS=16
ENV RENV_PATHS_CACHE=/renv_cache
ENV RENV_DOWNLOAD_FILE_METHOD=libcurl
ENV RENV_CONFIG_PAK_ENABLED=FALSE

# Ensure V8 can always find its library
ENV DOWNLOAD_STATIC_LIBV8=1

# ----------------------------
# Install system dependencies
# ----------------------------
# Merged dependencies from project and user reference
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    openssh-client \
    curl \
    ca-certificates \
    tar \
    gzip \
    libv8-dev \
    libxt6 \
    libgl1 \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    libproj-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libspatialite-dev \
    libcairo2-dev \
    libxt-dev \
    libmagick++-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    libpng-dev \
    libjpeg-dev \
    libicu-dev \
    libabsl-dev \
    libfftw3-dev \
    libnlopt-dev \
    libgit2-dev \
    tree \
    gdal-bin \
    fonts-dejavu \
    fontconfig \
    default-jdk \
    r-cran-rjava \
    && R CMD javareconf \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /renv_cache

# ----------------------------
# Copy renv lockfile
# ----------------------------
COPY renv.lock renv.lock

# ----------------------------
# Install renv and restore packages
# ----------------------------
RUN echo 'options(repos = c(CRAN="https://packagemanager.posit.co/cran/latest"))' \
    >> /usr/local/lib/R/etc/Rprofile.site

# Install specific renv version
RUN Rscript -e 'install.packages("remotes")' \
 && Rscript -e 'remotes::install_version("renv", version = "1.2.3")'

RUN Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE, exclude = "renv")'
