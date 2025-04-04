# Use the official Rocker image for R 4.4.0
FROM rocker/r-ver:4.4.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    curl \
    libblas-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libicu-dev \
    liblapack-dev \
    liblzma-dev \
    libfribidi-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Copy your repository files into the container
COPY . /app

# Install renv and restore the R package library from renv.lock
RUN Rscript -e "install.packages('renv', repos='https://cloud.r-project.org')" && \
    Rscript -e "renv::restore(prompt = FALSE)"

# Install R package dependencies
RUN Rscript -e "install.packages('devtools', repos='https://cloud.r-project.org')" && \
    Rscript -e "devtools::install_github('YuLab-SMU/ggtree@c17773c973d6c4036ee3af40a3957fb74d8ee9ff')" && \
    Rscript -e 'devtools::install_github("mmore500/mlscluster@a8e14c19c1ac6a75c539b898bf1c83216b613b1f")'

# Run the R script and consolidate outputs
CMD ["bash", "-c", "exec Rscript R/app.R \"$@\"", "--"]
