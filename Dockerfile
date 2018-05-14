FROM continuumio/anaconda

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  libgsl-dev \
  libncurses-dev \
  libz-dev

RUN conda create -n partis python=2.7 anaconda
RUN /bin/bash -c "source activate partis && conda install -y biopython pandas psutil pysam scons seaborn zlib"
RUN /bin/bash -c "source activate partis && conda install -y -c biocore mafft"
RUN /bin/bash -c "source activate partis && pip install colored-traceback dendropy==3.12.3"
RUN git clone https://github.com/psathyrella/partis.git
WORKDIR /partis
RUN /bin/bash -c "source activate partis && ./bin/build.sh"
RUN /bin/bash -c "source activate partis && conda install -y r-essentials"
RUN /bin/bash -c "source activate partis && unset R_LIBS_SITE && R --vanilla --slave -e 'install.packages(\"TreeSim\", repos=\"http://cran.rstudio.com/\")'"
