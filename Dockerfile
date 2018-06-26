FROM continuumio/anaconda

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
RUN conda install -y python=2.7
RUN conda install -y biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn
RUN conda install -y -c biocore mafft
RUN pip install colored-traceback dendropy==3.12.3
COPY . /partis
WORKDIR /partis
RUN ./bin/build.sh
CMD ./test/test.py --quick
