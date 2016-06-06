FROM matsengrp/cpp

# ----------------------------------------------------------------------------------------
# java stuff
RUN grep '8\.[0-9]' /etc/debian_version  # if it's not 8.x (jessie) we need to change the 'trusty' below
RUN echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" >> /etc/apt/sources.list
RUN echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" >> /etc/apt/sources.list
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886
RUN apt-get update #2>&1 | tee /tmp/apt.err && ! grep -q '^[WE]' /tmp/apt.err  # shenanigans to induce failure since apt-get update is unduly tolerant of failures
RUN echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
# ----------------------------------------------------------------------------------------
RUN apt-get install -y \
    astyle \
    oracle-java7-installer \
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
