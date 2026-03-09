# aangeloo/emu 
FROM continuumio/miniconda3:latest

#
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    tar \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y emu osfclient && \
    conda clean -afy

# Download the Emu databases
ENV EMU_DATABASE_DIR=/opt/emu_database
RUN mkdir -p ${EMU_DATABASE_DIR} && \
    cd ${EMU_DATABASE_DIR} && \
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar && \
    mkdir emu && tar -xvf emu.tar -C emu && rm emu.tar && \
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/rdp.tar && \
    mkdir rdp && tar -xvf rdp.tar -C rdp && rm rdp.tar && \
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/unite-all.tar && \
    mkdir unite-all && tar -xvf unite-all.tar -C unite-all && rm unite-all.tar


