FROM rocker/rstudio:4.6

# Script for to set up multiple users
# COPY batch_user_creation.sh /rocker_scripts/experimental/batch_user_creation.sh
# CMD [ "sh", "-c", "/rocker_scripts/experimental/batch_user_creation.sh && /init" ]

# Install system dependencies required for R packages
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    openssh-client \
    curl \
    ca-certificates \
    tar \
    gzip \
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
    libhdf5-dev \
    libnode-dev \
    tree \
    cmake \
    gdal-bin \
    fonts-dejavu \
    fontconfig \
    default-jdk \
    r-cran-rjava \
    && R CMD javareconf \
    && rm -rf /var/lib/apt/lists/*

COPY renv.lock ./

ENV RENV_INSTALL_THREADS=16
ENV RENV_DOWNLOAD_FILE_METHOD=libcurl
ENV RENV_CONFIG_PAK_ENABLED=FALSE

RUN Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
RUN Rscript -e 'remotes::install_version("renv", version = "1.2.3", repos="https://cloud.r-project.org")'

# install xgboost first (required by BayesSpace)
RUN Rscript -e 'install.packages("xgboost", repos = c("https://dmlc.r-universe.dev", "https://cloud.r-project.org"))'

# restore everything except renv, xgboost, and the two that collide in staging
RUN Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE, exclude = c("renv", "xgboost", "BPCells", "monocle3"))'

# install the stragglers explicitly (BPCells must build before monocle3)
RUN rm -rf /tmp/Rtmp*/renv-graph-staging-* 2>/dev/null; \
    Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE, packages = "BPCells")'
RUN Rscript -e 'renv::restore(lockfile="renv.lock", prompt=FALSE, packages = "monocle3")'

EXPOSE 8787
CMD ["/init"]