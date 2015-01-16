FROM matsengrp/cpp

# make ham
COPY . /ham
RUN cd /ham && \
    scons test
