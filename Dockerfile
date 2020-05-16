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
  libz-dev

RUN conda update -y conda
RUN conda install -y -cbioconda -cbiocore python=2.7 biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn mafft # -cbioconda is for pysam, -cbiocore is for mafft
RUN pip install colored-traceback dendropy==4.4.0

# ----------------------------------------------------------------------------------------
# this block is only necessary to install R for simulation (where we use the tree sim packages to generate initial trees -- alternatively you can use --input-simulation-treefname and skip all of this), and is necessary because if you install r with conda it replaces your entire compiler toolchain (e.g. screws up where your system looks for gcc) so you can't really compile anything afterwards
RUN apt-get install -y dirmngr apt-transport-https ca-certificates software-properties-common gnupg2
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/'  # this will have to change when anaconda changes debian releases, but the tree sim packages require a more recent version of r than is in the debian repos
RUN apt-get update && apt-get install -y r-base
RUN R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM"), repos="http://cran.rstudio.com/")'

# ----------------------------------------------------------------------------------------
COPY . /partis
WORKDIR /partis
RUN ./bin/build.sh with-simulation  # run with no arguments if you don't care about simulation, and it'll skip a compilation step (which is SLOOOOOOOOOOOW)
CMD ./test/test.py --quick
