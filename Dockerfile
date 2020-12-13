FROM continuumio/anaconda:5.3.0

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  libboost-dev \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libyaml-cpp-dev \
  libyaml-dev \
  libz-dev \
  python-pyqt5

RUN conda install -y -cbioconda -cbiocore python=2.7 biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn mafft # -cbioconda is for pysam, -cbiocore is for mafft
RUN conda update -y numpy  # the previous command downgrades numpy (I'm not sure why), which breaks the seaborn import
RUN pip install colored-traceback dendropy==4.4.0
COPY . /partis
WORKDIR /partis
RUN ./bin/build.sh
CMD ./test/test.py --quick
