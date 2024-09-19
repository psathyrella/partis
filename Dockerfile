FROM mambaorg/micromamba:1.5.6
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root
RUN apt-get --allow-releaseinfo-change update && apt-get install -y \
  automake \
  build-essential \
  cmake \
  emacs \
  fonts-lato \
  git \
  less \
  libboost-dev \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libyaml-cpp-dev \
  libyaml-dev \
  libz-dev \
  ncbi-blast+ \
  python3-pyqt5 \
  vim \
  xvfb
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -c conda-forge -c bioconda -c biocore -c etetoolkit \
    python biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn \
    mafft circlify ete3
RUN pip install colored-traceback dendropy joypy levenshtein
ENV XDG_RUNTIME_DIR=/tmp/xdgrd
RUN mkdir $XDG_RUNTIME_DIR
RUN chmod 0700 $XDG_RUNTIME_DIR
COPY --chown=$MAMBA_USER:$MAMBA_USER_GID . /partis
WORKDIR /partis
ENV TERM linux
ENV TERMINFO /etc/terminfo
RUN ./bin/build.sh
ENV PATH /opt/conda/bin:/opt/conda/condabin:$PATH  # necessary for singularity
CMD ./test/test.py --quick
