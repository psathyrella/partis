FROM continuumio/anaconda:5.3.0

RUN apt-get update && apt-get install -y \
  build-essential \
  libboost-dev \
  libgsl-dev \
  libncurses-dev \
  libxt6 \
  libyaml-cpp-dev \
  libyaml-dev \
  libz-dev

RUN apt-get install -y \
  cmake  # only for simulation

RUN conda update -y conda
RUN conda install -y -cbioconda -cbiocore python=2.7 biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn mafft # -cbioconda is for pysam, -cbiocore is for mafft
RUN conda install -y cmake r-essentials  # only for simulation
RUN unset R_LIBS_SITE  # only for simulation
RUN R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM"), repos="http://cran.rstudio.com/")'  # only for simulation
RUN pip install colored-traceback dendropy==4.4.0
COPY . /partis
WORKDIR /partis
RUN ./bin/build.sh
CMD ./test/test.py --quick
