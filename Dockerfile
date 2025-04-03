# Use the official Rocker image for R 3.6.3
FROM rocker/r-ver:3.6.3

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

# Install R package dependencies
RUN ./install_packages.sh

# Run the R script and consolidate outputs
CMD ["bash", "-c", "exec Rscript R/app.R \"$@\"", "--"]
