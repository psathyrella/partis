FROM matsengrp/cpp

# ----------------------------------------------------------------------------------------
RUN apt-get install -y \
    astyle \
    libgsl0ldbl \
    libgsl0-dev \
    libncurses5-dev \
    libxml2-dev \
    libxslt1-dev \
    python-scipy \
    r-base \
    zlib1g-dev
RUN pip install \
    numpy \
    scipy \
    matplotlib \
    pandas \
    biopython \
    cython \
    dendropy==3.12.3 \
    pysam \
    pyyaml \
    seaborn
RUN R --vanilla --slave -e 'install.packages("TreeSim", repos="http://cran.rstudio.com/")'

COPY . /partis
WORKDIR /partis
CMD ./bin/build.sh && ./test/test.py --quick
