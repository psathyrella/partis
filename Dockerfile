FROM matsengrp/cpp

# make ham
COPY . /ham
CMD cd /ham && \
    scons test
