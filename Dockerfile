FROM matsengrp/cpp

# Java bit copied from https://github.com/jplock/docker-oracle-java7
RUN sed 's/main$/main universe/' -i /etc/apt/sources.list
RUN apt-get update && apt-get install -y software-properties-common python-software-properties
RUN add-apt-repository ppa:webupd8team/java -y
RUN apt-get update
RUN echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
RUN apt-get install -y \
    oracle-java7-installer \
    libncurses5-dev \
    libroot-bindings-python-dev \
    libroot-graf2d-postscript5.34 \
    libxml2-dev \
    libxslt1-dev \
    python-scipy \
    zlib1g-dev
RUN pip install \
    beautifulsoup4 \
    biopython \
    cython \
    decorator \
    lxml \
    networkx \
    pysam \
    pyyaml

# set up auth
RUN mkdir -p /root/.ssh && \
    chmod 700 /root/.ssh
ADD bunnyhutch_id_rsa /root/.ssh/id_rsa
RUN chmod 600 /root/.ssh/id_rsa && \
    ssh-keyscan github.com >> /root/.ssh/known_hosts

# make
RUN git clone git@github.com:psathyrella/partis.git
WORKDIR /data/partis/packages/samtools/
RUN make && \
    mkdir /data/partis/_bin && \
    ln -s $PWD/samtools /data/partis/_bin
WORKDIR /data/partis/packages/ighutil/
RUN make -C clj && \
    pip install --user ./python
WORKDIR /data/partis/packages/ham/
RUN scons bcrham
WORKDIR /data/partis/
RUN mkdir -p /true/plots
RUN export PATH=/data/partis/_bin:$PATH && scons test
