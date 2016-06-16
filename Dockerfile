FROM matsengrp/cpp

# ----------------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    libgsl0ldbl \
    libgsl0-dev \
    libncurses5-dev \
    libxml2-dev \
    libxslt1-dev \
    r-base
RUN pip install numpy  # putting them on different lines allows docker's caching to defeat pip's slowness
RUN pip install scipy
RUN pip install cython
RUN pip install matplotlib
RUN pip install pandas
RUN pip install biopython
RUN pip install dendropy==3.12.3
RUN pip install pysam
RUN pip install pyyaml
RUN pip install seaborn

RUN R --vanilla --slave -e 'install.packages("TreeSim", repos="http://cran.rstudio.com/")'

COPY . /partis
WORKDIR /partis
CMD ./bin/build.sh && ./test/test.py --quick
