FROM matsengrp/cpp

# make ham
CMD git clone https://github.com/psathyrella/ham.git && \
    cd ham && \
    scons test
