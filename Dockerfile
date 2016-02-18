FROM matsengrp/cpp

# Java bit copied from https://github.com/jplock/docker-oracle-java7
RUN sed 's/main$/main universe/' -i /etc/apt/sources.list
RUN apt-get update && apt-get install -y software-properties-common python-software-properties
RUN add-apt-repository ppa:webupd8team/java -y
RUN sudo apt-get update 2>&1 | tee /tmp/apt.err && ! grep -q '^[WE]' /tmp/apt.err  # shenanigans to induce failure since apt-get update is unduly tolerant of failures
RUN echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
RUN apt-get install -y \
    astyle \
    oracle-java7-installer \
    libgsl0ldbl \
    libgsl0-dev \
    libncurses5-dev \
    libxml2-dev \
    libxslt1-dev \
    python-scipy \
    python-sklearn \
    r-base \
    zlib1g-dev
RUN pip install \
    beautifulsoup4 \
    biopython \
    cython \
    decorator \
    dendropy==3.12.3 \
    lxml \
    networkx \
    pysam \
    pyyaml \
    seaborn
RUN R --vanilla --slave -e 'install.packages("TreeSim", repos="http://cran.rstudio.com/")'


COPY . /partis
WORKDIR /partis
CMD ./bin/build.sh && export $PWD/packages/samtools && ./test/test.py --quick
