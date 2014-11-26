# make-ham
#
# VERSION       1.0
FROM dockerfile/python
MAINTAINER Frederick A. Matsen, matsen@fhcrc.org

# prepare ubuntu
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
    build-essential \
    git \
    libboost-dev \
    scons

# make ham
RUN git clone https://github.com/psathyrella/ham.git
RUN cd ham && scons
