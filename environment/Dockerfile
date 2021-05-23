# hash:sha256:feba7b537b78d7787c7fc1e42cd89fc1970cdaff25f00d76c180c147bfd0cfdb
FROM registry.codeocean.com/codeocean/r-base:3.6.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("here", "latest")' \
    && Rscript -e 'remotes::install_version("scales", "latest")' \
    && Rscript -e 'remotes::install_version("tidyverse", "latest")' \
    && Rscript -e 'remotes::install_version("viridis", "latest")'
