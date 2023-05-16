FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="fff07e47b375000d28fdaaf1a6523a4574cc710d780521a23f246442e164663a"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bedtools.yaml
#   prefix: /conda-envs/8bdc069892ecc8d0829bb8fd2bf472f0
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - bedtools =2.30
#     - gawk
#     - coreutils
#     - gzip
#     - grep
RUN mkdir -p /conda-envs/8bdc069892ecc8d0829bb8fd2bf472f0
COPY workflow/envs/bedtools.yaml /conda-envs/8bdc069892ecc8d0829bb8fd2bf472f0/environment.yaml

# Conda environment:
#   source: workflow/envs/damageprofiler.yaml
#   prefix: /conda-envs/6bba40e406a50c6ffbfc392272fa54e7
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - damageprofiler =1.1
RUN mkdir -p /conda-envs/6bba40e406a50c6ffbfc392272fa54e7
COPY workflow/envs/damageprofiler.yaml /conda-envs/6bba40e406a50c6ffbfc392272fa54e7/environment.yaml

# Conda environment:
#   source: workflow/envs/dedup.yaml
#   prefix: /conda-envs/170102d52a3e443be78cd15da657c0c8
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - dedup =0.12.8
#     - samtools =1.15
RUN mkdir -p /conda-envs/170102d52a3e443be78cd15da657c0c8
COPY workflow/envs/dedup.yaml /conda-envs/170102d52a3e443be78cd15da657c0c8/environment.yaml

# Conda environment:
#   source: workflow/envs/fastp.yaml
#   prefix: /conda-envs/fcfc12234984164b9fa400ca46295111
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - fastp =0.23.2
RUN mkdir -p /conda-envs/fcfc12234984164b9fa400ca46295111
COPY workflow/envs/fastp.yaml /conda-envs/fcfc12234984164b9fa400ca46295111/environment.yaml

# Conda environment:
#   source: workflow/envs/fgbio.yaml
#   prefix: /conda-envs/4ba0a61add4621c9de984d408692418f
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - fgbio =2.0.2
#     - samtools =1.15
RUN mkdir -p /conda-envs/4ba0a61add4621c9de984d408692418f
COPY workflow/envs/fgbio.yaml /conda-envs/4ba0a61add4621c9de984d408692418f/environment.yaml

# Conda environment:
#   source: workflow/envs/gatk.yaml
#   prefix: /conda-envs/4507de88da500c9fd301c961685808f3
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - gatk =3.8
RUN mkdir -p /conda-envs/4507de88da500c9fd301c961685808f3
COPY workflow/envs/gatk.yaml /conda-envs/4507de88da500c9fd301c961685808f3/environment.yaml

# Conda environment:
#   source: workflow/envs/genmap.yaml
#   prefix: /conda-envs/5d5997989b5c34b7766ab2d89b584bcb
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   
#   dependencies:
#     - genmap =1.3
RUN mkdir -p /conda-envs/5d5997989b5c34b7766ab2d89b584bcb
COPY workflow/envs/genmap.yaml /conda-envs/5d5997989b5c34b7766ab2d89b584bcb/environment.yaml

# Conda environment:
#   source: workflow/envs/mapping.yaml
#   prefix: /conda-envs/00a1fbfd086699548a489853868e06db
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   
#   dependencies:
#     - bwa =0.7.17
#     - samtools = 1.15
RUN mkdir -p /conda-envs/00a1fbfd086699548a489853868e06db
COPY workflow/envs/mapping.yaml /conda-envs/00a1fbfd086699548a489853868e06db/environment.yaml

# Conda environment:
#   source: workflow/envs/picard.yaml
#   prefix: /conda-envs/e757785681e1566fcaab0cfdbd29e30e
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - picard =2.27
RUN mkdir -p /conda-envs/e757785681e1566fcaab0cfdbd29e30e
COPY workflow/envs/picard.yaml /conda-envs/e757785681e1566fcaab0cfdbd29e30e/environment.yaml

# Conda environment:
#   source: workflow/envs/pruning.yaml
#   prefix: /conda-envs/c11e203b08591074197ec814f8e71d72
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python =3.11.*
#     - graph-tool =2.55
#     - pandas =2.0.1
RUN mkdir -p /conda-envs/c11e203b08591074197ec814f8e71d72
COPY workflow/envs/pruning.yaml /conda-envs/c11e203b08591074197ec814f8e71d72/environment.yaml

# Conda environment:
#   source: workflow/envs/python.yaml
#   prefix: /conda-envs/b8e21f01566bad016b4ddaf9230c096e
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   
#   dependencies:
#     - python =3.11.*
#     - numpy
RUN mkdir -p /conda-envs/b8e21f01566bad016b4ddaf9230c096e
COPY workflow/envs/python.yaml /conda-envs/b8e21f01566bad016b4ddaf9230c096e/environment.yaml

# Conda environment:
#   source: workflow/envs/qualimap.yaml
#   prefix: /conda-envs/05a66c43aafdbc848050b4a40694b344
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - qualimap =2.2.2a
RUN mkdir -p /conda-envs/05a66c43aafdbc848050b4a40694b344
COPY workflow/envs/qualimap.yaml /conda-envs/05a66c43aafdbc848050b4a40694b344/environment.yaml

# Conda environment:
#   source: workflow/envs/r-rectable.yaml
#   prefix: /conda-envs/9b95849aafbb6b9a45b5cd80ab66d516
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base =4.2.*
#     - r-data.table =1.14.8
#     - r-reactable =0.4.4
#     - r-rmarkdown =2.21
RUN mkdir -p /conda-envs/9b95849aafbb6b9a45b5cd80ab66d516
COPY workflow/envs/r-rectable.yaml /conda-envs/9b95849aafbb6b9a45b5cd80ab66d516/environment.yaml

# Conda environment:
#   source: workflow/envs/r.yaml
#   prefix: /conda-envs/2936ebf872f2a69a730831ab9786931c
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - r-base =4.2.*
#     - r-essentials =4.2
#     - r-svglite =2.1.1
#     - r-hmisc =5.1_0
RUN mkdir -p /conda-envs/2936ebf872f2a69a730831ab9786931c
COPY workflow/envs/r.yaml /conda-envs/2936ebf872f2a69a730831ab9786931c/environment.yaml

# Conda environment:
#   source: workflow/envs/repeatmasker.yaml
#   prefix: /conda-envs/3f6779567d73a0f63bb795e68a4d5a96
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - repeatmasker =4.1.2
#     - repeatmodeler =2.0.3
RUN mkdir -p /conda-envs/3f6779567d73a0f63bb795e68a4d5a96
COPY workflow/envs/repeatmasker.yaml /conda-envs/3f6779567d73a0f63bb795e68a4d5a96/environment.yaml

# Conda environment:
#   source: workflow/envs/samtools.yaml
#   prefix: /conda-envs/75d1a2e4dce42f48816190621a908088
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - samtools =1.15
RUN mkdir -p /conda-envs/75d1a2e4dce42f48816190621a908088
COPY workflow/envs/samtools.yaml /conda-envs/75d1a2e4dce42f48816190621a908088/environment.yaml

# Conda environment:
#   source: workflow/envs/shell.yaml
#   prefix: /conda-envs/d941761320cf9c9de9a22b634ba0383d
#   channels:
#     - conda-forge
#     - bioconda
#   
#   dependencies:
#     - gawk
#     - coreutils
#     - bc
#     - gzip
#     - grep
#     - perl
RUN mkdir -p /conda-envs/d941761320cf9c9de9a22b634ba0383d
COPY workflow/envs/shell.yaml /conda-envs/d941761320cf9c9de9a22b634ba0383d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.17.2/bio/picard/markduplicates/environment.yaml
#   prefix: /conda-envs/0e437010819dd7fdd387a5cdc1a87e5b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - picard =2.27
#     - samtools =1.15
#     - snakemake-wrapper-utils ==0.1.3
RUN mkdir -p /conda-envs/0e437010819dd7fdd387a5cdc1a87e5b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.17.2/bio/picard/markduplicates/environment.yaml /conda-envs/0e437010819dd7fdd387a5cdc1a87e5b/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/8bdc069892ecc8d0829bb8fd2bf472f0 --file /conda-envs/8bdc069892ecc8d0829bb8fd2bf472f0/environment.yaml && \
    mamba env create --prefix /conda-envs/6bba40e406a50c6ffbfc392272fa54e7 --file /conda-envs/6bba40e406a50c6ffbfc392272fa54e7/environment.yaml && \
    mamba env create --prefix /conda-envs/170102d52a3e443be78cd15da657c0c8 --file /conda-envs/170102d52a3e443be78cd15da657c0c8/environment.yaml && \
    mamba env create --prefix /conda-envs/fcfc12234984164b9fa400ca46295111 --file /conda-envs/fcfc12234984164b9fa400ca46295111/environment.yaml && \
    mamba env create --prefix /conda-envs/4ba0a61add4621c9de984d408692418f --file /conda-envs/4ba0a61add4621c9de984d408692418f/environment.yaml && \
    mamba env create --prefix /conda-envs/4507de88da500c9fd301c961685808f3 --file /conda-envs/4507de88da500c9fd301c961685808f3/environment.yaml && \
    mamba env create --prefix /conda-envs/5d5997989b5c34b7766ab2d89b584bcb --file /conda-envs/5d5997989b5c34b7766ab2d89b584bcb/environment.yaml && \
    mamba env create --prefix /conda-envs/00a1fbfd086699548a489853868e06db --file /conda-envs/00a1fbfd086699548a489853868e06db/environment.yaml && \
    mamba env create --prefix /conda-envs/e757785681e1566fcaab0cfdbd29e30e --file /conda-envs/e757785681e1566fcaab0cfdbd29e30e/environment.yaml && \
    mamba env create --prefix /conda-envs/c11e203b08591074197ec814f8e71d72 --file /conda-envs/c11e203b08591074197ec814f8e71d72/environment.yaml && \
    mamba env create --prefix /conda-envs/b8e21f01566bad016b4ddaf9230c096e --file /conda-envs/b8e21f01566bad016b4ddaf9230c096e/environment.yaml && \
    mamba env create --prefix /conda-envs/05a66c43aafdbc848050b4a40694b344 --file /conda-envs/05a66c43aafdbc848050b4a40694b344/environment.yaml && \
    mamba env create --prefix /conda-envs/9b95849aafbb6b9a45b5cd80ab66d516 --file /conda-envs/9b95849aafbb6b9a45b5cd80ab66d516/environment.yaml && \
    mamba env create --prefix /conda-envs/2936ebf872f2a69a730831ab9786931c --file /conda-envs/2936ebf872f2a69a730831ab9786931c/environment.yaml && \
    mamba env create --prefix /conda-envs/3f6779567d73a0f63bb795e68a4d5a96 --file /conda-envs/3f6779567d73a0f63bb795e68a4d5a96/environment.yaml && \
    mamba env create --prefix /conda-envs/75d1a2e4dce42f48816190621a908088 --file /conda-envs/75d1a2e4dce42f48816190621a908088/environment.yaml && \
    mamba env create --prefix /conda-envs/d941761320cf9c9de9a22b634ba0383d --file /conda-envs/d941761320cf9c9de9a22b634ba0383d/environment.yaml && \
    mamba env create --prefix /conda-envs/0e437010819dd7fdd387a5cdc1a87e5b --file /conda-envs/0e437010819dd7fdd387a5cdc1a87e5b/environment.yaml && \
    mamba clean --all -y
