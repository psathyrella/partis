FROM matsengrp/cpp

# Copied from https://github.com/jplock/docker-oracle-java7
RUN sed 's/main$/main universe/' -i /etc/apt/sources.list
RUN apt-get update && apt-get install -y software-properties-common python-software-properties
RUN add-apt-repository ppa:webupd8team/java -y
RUN apt-get update
RUN echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
RUN apt-get install -y \
    oracle-java7-installer \
    libxml2-dev \
    libxslt1-dev

RUN pip install \
    pysam \
    pyyaml \
    cython \
    networkx \
    decorator \
    lxml \
    bs4
RUN pip install --user ./python

# make ham
git clone git@github.com:psathyrella/partis.git
WORKDIR /data/partis/ham
RUN scons
WORKDIR /data/partis/samtools/
RUN make && ln -s $PWD/samtools ~/bin/
WORKDIR /data/partis/ighutil/
RUN make -C clj
WORKDIR /data/partis/
