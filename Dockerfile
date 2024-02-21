FROM mambaorg/micromamba:1.5.6
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root
RUN apt-get --allow-releaseinfo-change update && apt-get install -y \
  build-essential \
  cmake \
  emacs \
  libboost-dev \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libyaml-cpp-dev \
  libyaml-dev \
  libz-dev \
  ncbi-blast+ \
  python3-pyqt5 \
  vim
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -cconda-forge -cbioconda -cbiocore python biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn mafft # -cbioconda is for pysam, -cbiocore is for mafft
RUN pip install colored-traceback dendropy levenshtein
COPY --chown=$MAMBA_USER:$MAMBA_USER_GID . /partis
WORKDIR /partis
RUN ./bin/build.sh
CMD ./test/test.py --quick
