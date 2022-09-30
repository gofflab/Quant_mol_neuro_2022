FROM archrockfish/rstudio:2022.09.27_485
# FROM rocker/rstudio:4.0.4

LABEL org.opencontainers.image.source=https://github.com/gofflab/Quant_mol_neuro_2022
LABEL org.opencontainers.image.licenses=MIT

# Install system dependencies
RUN apt-get update
RUN apt-get install -y wget curl vim git zsh tldr && tldr ls

# R dependencies
RUN apt-get install -y zlib1g-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libxml2-dev libglpk-dev
RUN R -e "install.packages(\"pak\", repos = sprintf(\"https://r-lib.github.io/p/pak/stable/%s/%s/%s\", .Platform\$pkgType, R.Version()\$os, R.Version()\$arch))"
RUN R -e "pak::pkg_install(c('rmarkdown', 'tidyverse', 'here', 'ggridges', 'scRNAseq', 'scater', 'scran', 'bluster'))"

ADD docker/environment.yml /root
RUN /opt/conda/bin/mamba env update -f /root/environment.yml && \
    /opt/conda/bin/mamba clean --all -f -y

RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh" && \
    bash "Mambaforge-$(uname)-$(uname -m).sh" -b -p /opt/conda && \
    rm "Mambaforge-$(uname)-$(uname -m).sh" && \
    chown -R rstudio:rstudio /opt/conda

RUN su - rstudio -c '/opt/conda/bin/mamba init --all'

# RUN su - rstudio -c 'sh -c "$(wget -O- https://raw.githubusercontent.com/chaichontat/zsh-in-docker/6b63e78045933134ef7423de0ef9dedf4acb43d4/zsh-in-docker.sh)" -- \
#     -p https://github.com/zsh-users/zsh-autosuggestions \
#     -p https://github.com/zsh-users/zsh-completions \
#     -p https://github.com/zsh-users/zsh-history-substring-search \
#     -p https://github.com/zsh-users/zsh-syntax-highlighting'

# # RStudio config
# ADD docker/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json
# RUN chown -R rstudio:rstudio /home/rstudio/.config/rstudio

# ADD docker/.zshrc /home/rstudio/.zshrc
# RUN chown -R rstudio:rstudio /home/rstudio/.zshrc

# ADD docker/.p10k.zsh /home/rstudio/.p10k.zsh
# RUN chown -R rstudio:rstudio /home/rstudio/.p10k.zsh
