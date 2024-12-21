FROM mambaorg/micromamba:1.5.8 as app

ARG PIMA_VER="2.1.1"

# 'LABEL' instructions tag the image with metadata that might be important to the user
LABEL base.image="mambaorg/micromamba:1.5.8"
LABEL dockerfile.version="1"
LABEL software="pima"
LABEL software.version="${PIMA_VER}"
LABEL description="Plasmid, Integrations, Mutations, and Antibiotic resistance annotation pipeline"
LABEL maintainer="Will Overholt"
LABEL maintainer.email="woverholt@asrtinc.com"

USER root

## Need to build the file from the primary pima dir
# cd ....pima
# docker build -t local/pima:2.1.0 -f dockerbuild/Dockerfile .
# cd ....singularity_images
# apptainer build pima2.1.0.sif docker-daemon://local/pima:2.1.0

# include required pima files
ADD conda_recipe/environment.yml environment.yml
ADD Pima Pima
ADD setup.py setup.py
ADD README.md README.md

# build run environment
RUN apt-get update && apt-get install -y --no-install-recommends \
    locales \
    locales-all \
    libpango-1.0-0 libpangoft2-1.0-0 libharfbuzz-subset0 \
    wget \
    procps \
    ca-certificates && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN micromamba install --name base -c conda-forge -c bioconda -f environment.yml && \
    micromamba clean -a -f -y && \
    mkdir /data

ENV PATH="/opt/conda/bin/:$PATH" \
    LC_ALL=C.UTF-8

# install pima
RUN python -m pip install -vv --no-deps --ignore-installed --no-cache-dir . && \
    python -m pip install --no-deps --no-cache-dir dna_features_viewer si-prefix && \
    python -m pip install --force-reinstall --no-cache-dir --no-deps git+https://github.com/ponnhide/pyCircos.git

CMD pima --help

WORKDIR /data

FROM app as test

WORKDIR /test

# test installation
RUN pima --help && \
    pima --version

# prep pima
RUN pima --download