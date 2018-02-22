FROM matsengrp/cpp

# ----------------------------------------------------------------------------------------
RUN sed -i 's/httpredir/ftp.us/' /etc/apt/sources.list
# ape (required for TreeSim) doesn't work with the default debian r version; probably need to add the cran apt sources here
RUN apt-get update && apt-get install -y \
    libgsl0-dev \
    libncurses5-dev \
    libxml2-dev \
    libxslt1-dev \
    mafft \
    r-base
RUN pip install numpy  # putting them on different lines allows docker's caching to defeat pip's slowness
RUN pip install scipy
RUN pip install scikit-learn
RUN pip install matplotlib
RUN pip install pandas
RUN pip install biopython
RUN pip install dendropy==3.12.3
RUN pip install pysam
RUN pip install pyyaml
RUN pip install seaborn
RUN pip install colored_traceback
RUN pip install psutil

RUN R --vanilla --slave -e 'install.packages(c("TreeSim", "TreeSimGM", "bios2mds"), repos="http://cran.rstudio.com/")'

COPY . /partis
WORKDIR /partis
CMD ./bin/build.sh && ./test/test.py --quick
