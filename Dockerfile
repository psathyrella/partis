FROM dockerfile/python

# prepare ubuntu
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    libboost-dev \
    scons

# make ham
RUN git clone https://github.com/psathyrella/ham.git
WORKDIR /data/ham
RUN scons

# run tests
RUN scons test
