FROM condaforge/mambaforge:latest

LABEL maintainer="Vikas Peddu"
LABEL description="Docker image for long-read isoform analysis pipeline"

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH="/opt/conda/bin:$PATH"

# Ensure non-interactive apt installs and set default timezone
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        wget \
        bzip2 \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libncurses5-dev \
        samtools \
        bedtools \
        openjdk-11-jre-headless \
        r-base \
        r-base-dev \
        perl \
        vim \
        less \
        gawk \
        sed \
        grep \
        bc \
        time \
        python3-pip \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN mamba config --system --prepend channels conda-forge && \
    mamba config --system --prepend channels bioconda && \
    mamba config --system --set auto_update_conda false && \
    mamba config --system --set show_channel_urls true

RUN mamba create -y -n long_reads python=3.10 && \
    mamba create -y -n flair python=3.10 && \
    mamba create -y -n sqanti3 python=3.11 && \
    mamba create -y -n vep python=3.9

RUN /opt/conda/bin/mamba run -n long_reads pip install --upgrade pip && \
    /opt/conda/bin/mamba run -n long_reads pip install htseq pysam pandas numpy scipy biopython bcbio-gff pyyaml

RUN /opt/conda/bin/mamba install -y -n long_reads -c bioconda -c conda-forge minimap2 samtools bedtools seqtk pigz longshot subread bwa

RUN /opt/conda/bin/mamba run -n flair pip install --upgrade pip && \
    /opt/conda/bin/mamba install -y -n flair -c bioconda -c conda-forge flair bedtools samtools minimap2 biopython

RUN /opt/conda/bin/mamba run -n sqanti3 pip install --upgrade pip && \
    /opt/conda/bin/mamba install -y -n sqanti3 -c bioconda -c conda-forge sqanti3 cDNA_cupcake bedtools samtools minimap2 pandas numpy pysam biopython

RUN /opt/conda/bin/mamba install -y -n vep -c bioconda -c conda-forge ensembl-vep vcftools bcftools

RUN mkdir -p /opt/picard && cd /opt/picard && \
    curl -sSL -O "https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar" && \
    chmod +x picard.jar

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r)" > ~/.Rprofile && \
    Rscript -e "install.packages(c('tidyverse', 'ggplot2', 'data.table', 'devtools'), dependencies=TRUE, repos='http://cran.us.r-project.org')" || true

RUN mkdir -p /opt/bin && cat > /opt/bin/run_with_conda.sh << 'EOF'\n#!/bin/bash\nENV_NAME=$1\nshift\nsource /opt/conda/bin/activate $ENV_NAME\nexec "$@"\nEOF
RUN chmod +x /opt/bin/run_with_conda.sh

ENV PATH="/opt/bin:$PATH"

WORKDIR /workspace

CMD ["bash"]
