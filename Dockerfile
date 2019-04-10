FROM continuumio/anaconda

RUN cat /etc/os-release
RUN echo $PATH
RUN apt-get --help

RUN apt-get update && apt-get install -y \
  build-essential \
  libboost-dev \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libyaml-cpp-dev \
  libyaml-dev \
  libz-dev

RUN conda update -y conda
RUN conda install -y -cbioconda -cbiocore python=2.7 biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn mafft  # -cbioconda is for pysam, -cbiocore is for mafft
RUN pip install colored-traceback dendropy==4.0.0
COPY . /partis
WORKDIR /partis
RUN ./bin/build.sh
CMD ./test/test.py --quick
