FROM matsengrp/cpp

# Java bit copied from https://github.com/jplock/docker-oracle-java7
RUN sed 's/main$/main universe/' -i /etc/apt/sources.list
RUN apt-get update && apt-get install -y software-properties-common python-software-properties
RUN add-apt-repository ppa:webupd8team/java -y
RUN apt-get update
RUN echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
RUN apt-get install -y \
    oracle-java7-installer \
    libxml2-dev \
    zlib1g-dev \
    libxslt1-dev
RUN pip install \
    pysam \
    pyyaml \
    cython \
    networkx \
    decorator \
    lxml \
    beautifulsoup4

# set up
RUN mkdir -p /root/.ssh
ADD bunnyhutch_id_rsa /root/.ssh/id_rsa
RUN git config --global user.name bunnyhutch
RUN git config --global user.email ematsen+bunnyhutch@gmail.com
RUN ssh-keyscan github.com >> /root/.ssh/known_hosts

# make ham
RUN git clone https://github.com/matsengrp/partis.git
WORKDIR /data/partis/packages/ham/
RUN scons
WORKDIR /data/partis/packages/samtools/
RUN make && ln -s $PWD/samtools ~/bin/
WORKDIR /data/partis/packages/ighutil/
RUN make -C clj
WORKDIR /data/partis/
RUN pip install --user ./python
