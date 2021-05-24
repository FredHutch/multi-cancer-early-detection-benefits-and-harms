# hash:sha256:1d56c8640d047d8ab6defab3394dc9468851399ef8fee43dc4337631d9696128
FROM registry.codeocean.com/codeocean/r-base:3.6.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("scales", "1.1.1")' \
    && Rscript -e 'remotes::install_version("tidyverse", "1.3.0")' \
    && Rscript -e 'remotes::install_version("viridis", "0.5.1")'
