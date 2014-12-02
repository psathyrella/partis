FROM matsengrp/cpp

# make ham
RUN git clone https://github.com/psathyrella/ham.git
WORKDIR /data/ham
RUN scons

# run tests
RUN scons test
