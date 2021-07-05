FROM centos:7.6.1810
LABEL maintainer="Jumamurat Bayjanov"
LABEL version="0.1"

ENV MINICONDA_PATH /opt/miniconda3/bin/
ENV PATH="/opt/miniconda3/bin/:${PATH}"
ENV LC_ALL en_US.utf-8
ENV SHELL /bin/bash

RUN yum update -y && \
    yum install -y which libicu.x86_64 libicu-devel.x86_64 && \
    curl -o Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    chmod +x Miniconda3-latest-Linux-x86_64.sh && \
    ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    export PATH=$PATH:/opt/miniconda3/bin/ && \
    conda create --channel "conda-forge" --channel bioconda -n r_pkgs -y bioconductor-deseq2==1.24.0 r-lattice==0.20_38  r-pheatmap==1.0.12 bioconductor-edger==3.26.5 && \
    conda create --channel "conda-forge" --channel bioconda -n fastp -y fastp==0.20.0 && \
    conda create --channel "conda-forge" --channel bioconda -n multiqc -y multiqc==1.7 && \
    conda create --channel "conda-forge" --channel bioconda -n star -y star==2.7.1a && \
    conda create --channel "conda-forge" --channel bioconda -n rsem -y rsem==1.3.1 

RUN export PATH=$PATH:/opt/miniconda3/bin/ && \
    conda init -q bash && \
    bash -c "source ~/.bashrc" && \
    /bin/bash && \
    exec /bin/bash && \
    conda clean --tarballs && \
    rm ./Miniconda3-latest-Linux-x86_64.sh

#CMD bash -c "source ~/.bashrc" && exec "$@"

#    Author  JBayjan
#    Version v0.0.3
#    Purpose NGS pipeline
#
#%help
#    This image contains the following tools fastp==0.20.0, multiqc==1.7, star==2.7.1a, rsem==1.3.1.
#    And following R packages bioconductor-deseq2==1.24.0, r-lattice==0.20_38, r-pheatmap==1.0.12, bioconductor-edger==3.26.5.
#
# use this image by activating/deactivating corresponding environment for each tool.
#    In order to use tools activate & deactivate the corresponding environments of each tool:
#    E.g: You can put activate/deactivate lines at the beginning and at the end of you script for each step.
#    source activate rsem
#    ...do you work here...
#    conda deactivate
#    To use R packages activate and once finished deactivate the environment: r_pkgs.
#
